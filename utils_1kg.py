from inventory import *
from contextlib import contextmanager
from scipy.stats import norm, kstest, shapiro


def load_genome_names():
    """
    Loads a list of genome names from storage.
    """
    sample_names = load_sample_names()
    genome_names = [f"{name}_1" for name in sample_names] + [f"{name}_2" for name in sample_names]
    return genome_names


def load_sample_names():
    """
    Loads a list of sample names from storage.
    """
    dfs = load_df("european_samples")
    sample_names = dfs.name.tolist()
    return sample_names


def load_unrelated_sample_names():
    """
    Loads a list of unrelated European sample names from storage.
    """
    dfs = load_df("european_samples")
    unrelated_sample_names = dfs[~dfs.related].name.tolist()
    return unrelated_sample_names


def load_unrelated_genome_names():
    """
    Loads and formats genome names for unrelated European samples.
    
    This function:
    1. Loads a list of unrelated European sample names from storage
    2. Creates two genome names per sample (1 and 2) to represent diploid chromosomes
    
    Returns:
        list: A list of genome names in the format [name_1, name_2] for each European sample
    """
    # Load the base list of unrelated European sample names
    unrelated_sample_names = load_unrelated_sample_names()

    # Create two genome names per sample by appending _1 and _2 
    # This represents the two chromosomes in diploid cells
    unrelated_genome_names = [f"{name}_1" for name in unrelated_sample_names] + [f"{name}_2" for name in unrelated_sample_names]

    return unrelated_genome_names


def flip_zero_one(x):
    """
    Flips 0 to 1 and 1 to 0.
    Used to convert boolean of having alt allele to boolean of having ref allele.

    Args:
        x: Input value (0 or 1)

    Returns:
        The flipped value (1 if input is 0, 0 if input is 1)
    """
    return (x - 1) * -1


def flip_zero_two(x):
    """
    Flips 0 to 2 and 2 to 0.
    Used to convert boolean of having alt allele to boolean of having ref allele.

    Args:
        x: Input value (0 or 1)

    Returns:
        The flipped value (2 if input is 0, 0 if input is 2)
    """
    return (x - 2) * -1



def download_file(url, file_path):
    response = requests.get(url, stream=True)
    response.raise_for_status()
    with open(file_path, 'wb') as file:
        for chunk in response.iter_content(chunk_size=8192):
            file.write(chunk)


def check_md5(filepath, filepath_md5):
    """
    Check the MD5 checksum of a file.
    """
    # Calculate the MD5 checksum of the file
    with open(filepath, 'rb') as f:
        file_hash = hashlib.md5()
        while chunk := f.read(8192):
            file_hash.update(chunk)
    calculated_md5 = file_hash.hexdigest()

    # Read the expected MD5 checksum from the .md5 file
    with open(filepath_md5, 'r') as f:
        expected_md5 = f.read().strip().split()[0]

    # Compare the calculated and expected MD5 checksums
    if calculated_md5 != expected_md5:
        raise ValueError(f"MD5 checksum does not match for {filepath}. Expected {expected_md5}, got {calculated_md5}.")


def gunzip(filepath):
    """Decompress a .gz file."""
    output_filepath = filepath.removesuffix(".gz")
    with gzip.open(filepath, 'rb') as f_in:
        with open(output_filepath, 'wb') as f_out:
            f_out.write(f_in.read())


@contextmanager
def suppress_output():
    """Suppress stdout and stderr."""
    with open(os.devnull, 'w') as fnull:
        old_stdout = sys.stdout
        old_stderr = sys.stderr
        try:
            sys.stdout = fnull
            sys.stderr = fnull
            yield
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr

def curate_european_samples():
    """
    Filter samples from the 1000 Genomes Project data to select European samples,
    excluding Finnish populations and optionally removing children.
    
    The function:
    1. Loads pedigree and sample information
    2. Optionally filters out children based on parent IDs
    3. Joins pedigree and sample data
    4. Filters for European ancestry samples
    5. Excludes Finnish populations
    6. Saves the filtered sample list
    """
    # Load the unrelated and related genome index files
    dfu = pd.read_csv('data/unrelated.txt', sep='\t', skiprows=lambda x: x < 23)
    dfr = pd.read_csv('data/related.txt', sep='\t', skiprows=lambda x: x < 23)

    # Get the samples (for marking related vs unrelated)
    samples_unrelated = set(dfu.SAMPLE_NAME)  # The NYGS is stupid 
    samples_related = set(dfr.SAMPLE_ID)  # The NYGS is stupid 

    # I'm serious they are stupid
    samples_unrelated = {s.rstrip() for s in samples_unrelated}
    samples_related = {s.rstrip() for s in samples_related}

    # Load the samples file
    dfs = pd.read_csv("data/samples.txt", sep=" ")

    # Remove the person with African ancestry (HG01783) and their child
    dfs = dfs.loc[dfs.SampleID != "HG01783"]  # remove the person
    dfs = dfs.loc[dfs.FatherID != "HG01783"]  # remove the child 

    # "HG00144" was the mother of "HG00155", but "HG00144" is not in the samples file
    dfs.loc[dfs.SampleID == "HG00155", "MotherID"] = '0'

    # Select the European superpopulation
    dfs = dfs[dfs.Superpopulation == "EUR"]

    # Convert the sex column to strings
    dfs.Sex = dfs.Sex.map({1: 'male', 2: 'female'})

    # Mark related vs unrelated
    dfs["related"] = np.nan
    dfs.loc[dfs.SampleID.isin(samples_related), "related"] = 1
    dfs.loc[dfs.SampleID.isin(samples_unrelated), "related"] = 0

    assert dfs.related.isna().sum() == 0

    dfs.related = dfs.related.astype(bool)

    assert set(dfs.query("related").SampleID).issubset(samples_related)
    assert set(dfs.query("~related").SampleID).issubset(samples_unrelated)

    # Drop the superpopulation column
    assert is_all(dfs.Superpopulation, "EUR")
    dfs.drop(columns=["Superpopulation"], inplace=True)

    # Remove Finnish population
    dfs = dfs[dfs.Population != "FIN"]

    # Define column name mappings for cleaner, more consistent naming
    rename = {'SampleID': "name",
                'Sex': "sex", 
                'Population': "population",
                'FatherID': "father_name",
                'MotherID': "mother_name",
                'FamilyID': "family"}

    # Rename columns using the mapping
    dfs.rename(columns=rename, inplace=True)

    # Replace "0" values with None for parent fields to indicate no parent
    dfs.loc[dfs.father_name == "0", "father_name"] = None
    dfs.loc[dfs.mother_name == "0", "mother_name"] = None

    # Add boolean column indicating if sample is a child (has at least one parent)
    dfs["child"] = dfs.father_name.notna() | dfs.mother_name.notna()

    # Check for missing parent data in children
    # Count number of children with missing father/mother data
    n_deadbeat_dads, n_deadbeat_moms = dfs.query("child").father_name.isna().sum(), dfs.query("child").mother_name.isna().sum()

    # Verify that there are no children with missing parent data
    assert n_deadbeat_dads == n_deadbeat_moms == 0

    # Add a parent column
    parents = set(dfs.father_name).union(set(dfs.mother_name))
    dfs["parent"] = dfs.name.isin(parents)

    # Sort the data
    dfs.sort_values(["child", "population", "sex"], ascending=[False, True, False], inplace=True)
    dfs.reset_index(drop=True, inplace=True)

    trios = dfs.loc[dfs.child][["name", "father_name", "mother_name"]].copy()

    trios.reset_index(drop=True, inplace=True)

    name_to_child_name = defaultdict(lambda: None)
    name_to_trio = defaultdict(lambda: None)
    for i, row in trios.iterrows():
        # row["name"] is the name of the child
        name_to_trio[row["name"]] = i + 1
        name_to_trio[row.father_name] = i + 1
        name_to_trio[row.mother_name] = i + 1
        name_to_child_name[row.father_name] = row["name"]
        name_to_child_name[row.mother_name] = row["name"]

    dfs["trio"] = dfs.name.map(name_to_trio)
    # Convert trio dtype to Nullable Int64
    dfs.trio = dfs.trio.astype("Int64")

    dfs["child_name"] = dfs.name.map(name_to_child_name)

    dfs.sort_values(["trio", "child", "sex"], ascending=[True, False, False], inplace=True)

    dfs.reset_index(drop=True, inplace=True)

    # There are no 'related' samples that are not already in a trio
    assert dfs[(dfs.related) & (dfs.trio.isna())].empty
    # Just another way of writing the same thing
    assert dfs[(dfs.related) & (~dfs.parent) & (~dfs.child)].empty

    # There are no children that are also parents
    assert dfs[dfs.child & dfs.parent].empty

    # Reorder columns to put key identifying information first
    dfs = place_first(dfs, ["name", "sex", "population", "trio", "child", "parent", "child_name", "father_name", "mother_name", "family"])
    # Verify that all sample names are unique to avoid duplicates
    assert dfs.name.is_unique

    # Check that the missing mother is not in our data 
    assert "HG00144" not in dfs.name

    # Save the processed dataframe and sample names list to files
    # This allows the data to be used in other parts of the pipeline
    save_df(dfs, "european_samples")

    # Test the european samples
    check_european_samples()


def extract_european_genome(chrom):
    """
    Process a VCF file for specific samples and extract the European genomes, filter to their common variants. 
    
    Args:
        chrom (str): Chromosome number to process
    """
    # Define file paths using f-strings

    if chrom == "X":
        vcf_file = f"data/genomes/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    else:
        vcf_file = f"data/genomes/1kGP_high_coverage_Illumina.chr{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    output_dir = "data/genomes2"
    temp_file = f"{output_dir}/temp_chr{chrom}.vcf.gz"
    output_file = f"{output_dir}/processed_chr{chrom}.csv"
    
    # Create output directory if it doesn't exist
    ensure_dir(output_dir)
    
    # Read chosen samples
    european_sample_names = load_sample_names()
    
    # Create a string of samples for bcftools
    samples_str = ','.join(european_sample_names)
    
    # First, extract only the chosen samples and filter for AC_EUR_unrel > 0
    # The -i option filters for AC_EUR_unrel > 0
    # The -s option selects only the specified samples
    filter_command = f"""bcftools view -i 'AC_EUR_unrel > 0' \
        -s {samples_str} {vcf_file} -Oz -o {temp_file}"""
    
    # Execute the filtering command
    subprocess.run(filter_command, shell=True, check=True)
    
    # Now convert to CSV format with desired fields
    # Include CHROM, POS, ID, REF, ALT, and genotypes for all samples
    csv_command = f"""bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%AC_EUR_unrel\\t[%GT]\\n' \
        {temp_file} > {output_file}"""
    
    # Execute the CSV conversion command
    subprocess.run(csv_command, shell=True, check=True)
    
    # Clean up temporary file
    os.remove(temp_file)


def extract_genotype_matrix(chrom):
    """
    Extract genotype data from European samples. 
    
    This function reads a processed csv file for a given chromosome,
    extracts genotype data for all European samples, and 
    performs validation checks against the original VCF data.

    Args:
        chrom (str): Chromosome number/name to process (e.g. "21", "X")

    The function:
    1. Loads European sample names
    2. Reads and validates variant data (position, ref/alt alleles)
    3. Creates and validates CPRA identifiers
    4. Extracts individual genotypes for all samples
    5. Performs random validation checks against original VCF
    """

    # Load the list of European sample names 
    european_sample_names = load_sample_names()

    # Construct filepath to the processed VCF data for this chromosome
    filepath = f"data/genomes2/processed_chr{chrom}.csv"

    # Define column names for the VCF data
    # acu (ac_eur_unrel) represents allele count in unrelated European samples
    vcf_cols = ['chrom', 'pos', 'cpra', 'ref', 'alt', 'acu', 'genotypes']

    # Read the raw VCF data into a pandas DataFrame
    # The file is tab-separated and we specify the column names defined above
    df = pd.read_csv(filepath, sep='\t', names=vcf_cols)

    print(chrom, df.shape[0])

    # Verify that the CPRA (Chrom:Pos:Ref:Alt) identifiers are unique
    # This is important as CPRA serves as a unique identifier for variants
    assert df.cpra.is_unique

    # Create our own CPRA identifier by concatenating chromosome (without 'chr' prefix),
    # position, reference allele, and alternate allele
    df["my_cpra"] = df.chrom.map(lambda x: x[3:]) + ":" + df.pos.astype(str) + ":" + df.ref + ":" + df.alt

    # Compare our CPRA with the one provided in the VCF
    # This helps identify any discrepancies in variant representation
    df["cpra_match"] = df.my_cpra == df.cpra

    # I cannot verify that our CPRA identifiers are also unique (it wasn't for chromosome 1)
    # assert df.my_cpra.is_unique

    # Note on CPRA discrepancies:
    # Through testing, we found that our CPRA format matches gnomAD when differences occur
    # The main source of differences appears to be in position numbering
    # The exact reason for these positional differences remains unclear

    # Process the genotype data for all samples
    # The genotypes string contains concatenated genotypes for all samples
    # Each genotype takes 3 characters (e.g., "0|0" or "1|0")
    genotype_dict = {
        name: df.genotypes.map(lambda x: x[i*3:i*3+3])
        for i, name in enumerate(european_sample_names)
    }

    # Add all individual sample genotypes to the DataFrame in one operation
    # This is more efficient than adding columns one at a time
    df = pd.concat([df, pd.DataFrame(genotype_dict)], axis=1)
    
    # Test agreement between our data processed data and the original vcf.gz file
    # for _ in range(1_000):
    #     name = np.random.choice(european_sample_names)
    #     test_agreement_of_genotypes(df, chrom, name, min_acu=50)  # Only test variants with European allele count >= 50

    # Remove the raw genotypes string column since we've already extracted individual genotypes
    df.drop(columns=["genotypes"], inplace=True)

    # Verify chromosome consistency
    assert is_all(df.chrom, f"chr{chrom}")

    # Remove chromosome column since it's redundant (all rows are the same chromosome)
    df.drop(columns=['chrom'], inplace=True)

    # Verify no missing values in the dataset
    assert df.isna().sum().sum() == 0

    # Create dictionaries to store the first allele (before '|') for each sample
    # This represents the allele inherited from one parent
    phased_binary_1_dict = {
        f"{name}_1": df[name].map(lambda x: x[0]).astype(int)
        for name in european_sample_names
    }   

    # Create dictionaries to store the second allele (after '|') for each sample
    # This represents the allele inherited from the other parent
    phased_binary_2_dict = {
        f"{name}_2": df[name].map(lambda x: x[2]).astype(int)
        for name in european_sample_names
    }

    # Add the separated alleles as new columns to the dataframe
    df = pd.concat([df, pd.DataFrame(phased_binary_1_dict), pd.DataFrame(phased_binary_2_dict)], axis=1)

    # Verify that our allele separation was correct by comparing with original genotypes
    for name in european_sample_names:
        assert is_all(df[name].map(lambda x: x[0]).astype(int) == df[f"{name}_1"])
        assert is_all(df[name].map(lambda x: x[2]).astype(int) == df[f"{name}_2"])

    # Remove the original genotype columns since we now have them separated into alleles
    df.drop(columns=european_sample_names, inplace=True)

    # Names of the diploid genomes
    genome_names = load_genome_names()
    
    # Calculate raw allele count (acr) by summing across all genomes
    # This represents the total number of alternate alleles observed in our dataset
    df["acr"] = df[genome_names].sum(axis=1)

    # Reorder columns to group related fields together
    # First 5 columns (pos, cpra, ref, alt, acu) followed by acr, then remaining columns
    column_reorder = df.columns[:5].tolist() + ["acr"] + df.columns[5:-1].tolist()
    df = df[column_reorder]

    # Filter variants to keep only those present in our population
    # acr >= 1 means the variant was observed at least once
    df = df.loc[df.acr >= 1]

    # Save processed data to CSV file
    output_dir = "data/genomes3"
    ensure_dir(output_dir)
    output_file = f"{output_dir}/processed_chr{chrom}.csv"
    df.to_csv(output_file, index=False)  # Save without row indices


def extract_x_chr_genotype_matrix():
    chrom = "X"

    # Load the list of European sample names 
    european_sample_names = load_sample_names()

    # Construct filepath to the processed VCF data for this chromosome
    filepath = f"data/genomes2/processed_chr{chrom}.csv"

    # Define column names for the VCF data
    # acu (ac_eur_unrel) represents allele count in unrelated European samples
    vcf_cols = ['chrom', 'pos', 'cpra', 'ref', 'alt', 'acu', 'genotypes']

    # Read the raw VCF data into a pandas DataFrame
    # The file is tab-separated and we specify the column names defined above
    df = pd.read_csv(filepath, sep='\t', names=vcf_cols)

    # Verify that the CPRA (Chrom:Pos:Ref:Alt) identifiers are unique
    # This is important as CPRA serves as a unique identifier for variants
    if chrom != "X":
        assert df.cpra.is_unique

    # Create our own CPRA identifier by concatenating chromosome (without 'chr' prefix),
    # position, reference allele, and alternate allele
    df["my_cpra"] = df.chrom.map(lambda x: x[3:]) + ":" + df.pos.astype(str) + ":" + df.ref + ":" + df.alt

    # Compare our CPRA with the one provided in the VCF
    # This helps identify any discrepancies in variant representation
    df["cpra_match"] = df.my_cpra == df.cpra

    # I cannot verify that our CPRA identifiers are also unique (it wasn't for chromosome 1)
    # assert df.my_cpra.is_unique

    # Note on CPRA discrepancies:
    # Through testing, we found that our CPRA format matches gnomAD when differences occur
    # The main source of differences appears to be in position numbering
    # The exact reason for these positional differences remains unclear

    # Process the genotype data for all samples
    # The genotypes string contains concatenated genotypes for all samples
    # Each genotype takes 3 characters (e.g., "0|0" or "1|0")

    def process_x_genotype(genotypes):
        # Split genotypes string into individual genotypes
        genotypes_list = []
        i = 0
        while i < len(genotypes):
            # Check if next 3 chars contain '|' pattern
            if i+2 < len(genotypes) and genotypes[i+1] == '|':
                genotypes_list.append(genotypes[i:i+3])
                i += 3
            else:
                genotypes_list.append(genotypes[i])
                i += 1
        assert len(genotypes_list) == len(european_sample_names)
        return tuple(genotypes_list)

    # Process all genotypes at once with a single map operation
    all_genotypes = df.genotypes.map(process_x_genotype)

    # Create dictionary of genotypes for each sample by position
    genotype_dict = {
        name: [genotypes[i] for genotypes in all_genotypes]
        for i, name in enumerate(european_sample_names)
    }

    # Add all individual sample genotypes to the DataFrame in one operation
    df = pd.concat([df, pd.DataFrame(genotype_dict)], axis=1)

    # Test agreement between our data processed data and the original vcf.gz file
    # for _ in range(1_000):
    #     name = np.random.choice(european_sample_names)
    #     test_agreement_of_genotypes(df, chrom, name, min_acu=50)  # Only test variants with European allele count >= 50


    # Remove the raw genotypes string column since we've already extracted individual genotypes
    df.drop(columns=["genotypes"], inplace=True)

    # Verify chromosome consistency
    assert is_all(df.chrom, f"chr{chrom}")

    # Remove chromosome column since it's redundant (all rows are the same chromosome)
    df.drop(columns=['chrom'], inplace=True)

    # Verify no missing values in the dataset
    assert df.isna().sum().sum() == 0

    # Create dictionaries to store the first allele (before '|') for each sample
    # This represents the allele inherited from one parent
    phased_binary_1_dict = {
        f"{name}_1": df[name].map(lambda x: x[0]).astype(int)
        for name in european_sample_names
    }   

    # Note: for the phased PAR regions, we cannot definitively determine which phase represents which sex chromosome (X vs Y)
    # we need to email the authors for this information

    # We are implicity assuming here that the first allele is the X chromosome and the second allele is the Y chromosome for males

    # We impute second X chromosome sites in males with -1, we can change this later 

    # Create dictionaries to store the second allele (after '|') for each sample
    # This represents the allele inherited from the other parent
    phased_binary_2_dict = {
        f"{name}_2": df[name].map(lambda x: x[2] if len(x) > 1 else -1).astype(int)
        for name in european_sample_names
    }

    # Add the separated alleles as new columns to the dataframe
    df = pd.concat([df, pd.DataFrame(phased_binary_1_dict), pd.DataFrame(phased_binary_2_dict)], axis=1)

    # Verify that our allele separation was correct by comparing with original genotypes
    for name in european_sample_names:
        assert is_all(df[name].map(lambda x: x[0]).astype(int) == df[f"{name}_1"])
        assert is_all(df[name].map(lambda x: x[2] if len(x) > 1 else -1).astype(int) == df[f"{name}_2"])

    # Remove the original genotype columns since we now have them separated into alleles
    df.drop(columns=european_sample_names, inplace=True)

    # Names of the diploid genomes
    genome_names = load_genome_names()

    # Calculate raw allele count (acr) by summing across all genomes
    # This represents the total number of alternate alleles observed in our dataset
    df["acr"] = df[genome_names].replace(-1, 0).sum(axis=1)

    # Reorder columns to group related fields together
    # First 5 columns (pos, cpra, ref, alt, acu) followed by acr, then remaining columns
    column_reorder = df.columns[:5].tolist() + ["acr"] + df.columns[5:-1].tolist()
    df = df[column_reorder]

    # Filter variants to keep only those present in our population
    # acr >= 1 means the variant was observed at least once
    df = df.loc[df.acr >= 1]

    # Save processed data to CSV file
    output_dir = "data/genomes3"
    ensure_dir(output_dir)
    output_file = f"{output_dir}/processed_chr{chrom}.csv"
    df.to_csv(output_file, index=False)  # Save without row indices


def download_catalog_score(score_name):
    """
    Downloads a score file from the PGS Catalog, harmonized to GRCh38.
    Checks that the MD5 checksums match the downloaded files.
    """
    if os.path.exists(f"data/catalog_scores/{score_name}.txt.gz"):
        return
    time.sleep(4)  # Sleep to avoid rate limiting

    score_file_url = f"https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{score_name}/ScoringFiles/Harmonized/{score_name}_hmPOS_GRCh38.txt.gz"
    output_dir = "data/catalog_scores"
    ensure_dir(output_dir)
    score_filepath = f"{output_dir}/{score_name}.txt.gz"
    md5_file_url = score_file_url + ".md5"
    md5_filepath = score_filepath + ".md5"

    # Download the score and md5 files
    download_file(score_file_url, score_filepath)
    download_file(md5_file_url, md5_filepath)

    # Verify the checksum
    check_md5(score_filepath, md5_filepath)


def process_catalog_score(score_name):
    """
    Load and preprocess a polygenic score file from the PGS Catalog.

    This function loads a polygenic score file from the PGS Catalog, performs validation checks,
    and standardizes the data format. The score files contain SNP information and effect sizes
    used for calculating polygenic risk scores.

    Args:
        score (str): The PGS Catalog score ID (e.g. 'PGS004211')

    Returns:
        pd.DataFrame: A preprocessed DataFrame containing the score data with columns:
            - chrom: Chromosome number (str, e.g. '1', '2', etc.)
            - pos: Genomic position (int, e.g. 818725)
            - rsid: SNP rsID if available (e.g. 'rs12184325')
            - effect_allele: Effect allele (e.g. 'C', 'A', etc.)
            - other_allele: Other/reference allele (e.g. 'T', 'G', etc.)
            - effect_weight: Effect size/weight (float, e.g. 0.007617)
            - variant_description: Additional variant info if available
            - hm_source: Source of harmonization (e.g. 'ENSEMBL')

    Raises:
        ValueError: If the score file has no rows after preprocessing
        ValueError: If neither 'other_allele' nor 'hm_inferOtherAllele' columns are present

    Example output:
        effect_allele  other_allele  effect_weight  ...  rsid        chrom    pos
        C             T             0.007617       ...  rs12184325   1        818725
        A             G            -0.008576       ...  rs2710890    1        1023525
    """

    if os.path.exists(f"data/scores/{score_name}.csv"):
        return

    # Define file paths
    filepath = f"data/catalog_scores/{score_name}.txt.gz"
    filepath_decompressed = filepath.removesuffix(".gz")

    # Check if file needs to be decompressed
    if not os.path.isfile(filepath_decompressed):
        gunzip(filepath)  # Decompress file
    
    # Load the score file
    dft = pd.read_csv(filepath_decompressed, sep='\t', comment='#', low_memory=False)
    
    # Standardize column names
    dft = dft.rename({"hm_chr": "chrom", "hm_pos": "pos", "hm_rsID": "rsid"}, axis=1)
    
    # Remove redundant columns if they exist
    for col in ["rsID", "chr_name", "chr_position"]:
        if col in dft.columns:
            dft.drop(columns=[col], inplace=True)
    
    # Handle missing chromosome values
    num_nan = dft.chrom.isna().sum()
    if num_nan > 0:
        dft = dft.dropna(subset=["chrom"])
    
    # Convert chromosome to string type and validate
    if dft.chrom.dtype == float:
        dft.chrom = dft.chrom.astype(int).astype(str)
        
    dft.chrom = dft.chrom.astype(str)
    if not set(dft.chrom).issubset(CHROMOSOMES):
        dft = dft[dft.chrom.isin(CHROMOSOMES)]
    
    # Handle other allele information
    if "hm_inferOtherAllele" in dft.columns:
        if dft.hm_inferOtherAllele.isna().all():
            dft.drop(columns=["hm_inferOtherAllele"], inplace=True)
            
    if "other_allele" not in dft.columns:
        if "hm_inferOtherAllele" not in dft.columns:
            raise ValueError(f"The DataFrame must have an 'other_allele' or 'hm_inferOtherAllele' column, but has {dft.columns}.")
        dft["other_allele"] = dft.hm_inferOtherAllele
        
    if dft.other_allele.isna().sum() > 0:
        dft = dft.dropna(subset=["other_allele"])

    # Final validation check
    if dft.shape[0] == 0:
        raise ValueError(f"No rows in {score_name} score.")

    # Delete this possible columns (we can incorporate them later if we want)
    for col in ["variant_description", "hm_source"]:
        if col in dft.columns:
            dft.drop(columns=[col], inplace=True)
    output_dir = "data/scores"
    ensure_dir(output_dir)
    output_file = f"{output_dir}/{score_name}.csv"
    assert one_per(dft, ['chrom', 'pos'])
    dft.to_csv(output_file, index=False)


def intersect_score(score_name: str) -> None:
    """
    Intersect a polygenic score with the European samples' genotype matrices. 

    Run time: about 10-15 minutes depending on how big the score is.

    This function takes a polygenic score DataFrame and matches its variants to the samples' genotype matrices,
    handling various edge cases like missing alleles, ambiguous variants, and multiallelic sites.
    Results are saved to separate files based on variant categories.

    Args:
        dft (pd.DataFrame): Polygenic score DataFrame containing variants and effect sizes
        score (str): Score identifier used for output file naming

    The function processes variants in these categories:
    - Selected: Variants that pass all filters and will be used
    - Missing: Variants where effect allele doesn't match reference or alternate allele
    - Ambiguous: Variants with mismatched reference-alternate and effect-other allele pairs
                 (effect allele is represented but not the other allele)
    - Multiallelic ambiguous: Multiallelic variants that can't be uniquely resolved into a matching bi-allelic pair
    """
    if os.path.exists(f"data/score-genomes/{score_name}"):
        return
        
    dft = pd.read_csv(f"data/scores/{score_name}.csv")

    assert {'effect_allele', 'other_allele', 'chrom', 'pos', 'effect_weight'}.issubset(dft.columns)

    # Select the chromosomes to use in the polygenic score
    chromosomes_score = set(dft.chrom)
    chromosomes_score = sort_chrom(chromosomes_score)

    # Create lists to store DataFrames for each category of sites
    dts_list = []              # Sites that pass all filters and will be used
    missing_list = []          # Sites where effect allele doesn't match ref or alt (missing)
    ambiguous_list = []        # Sites with mismatched ref-alt and effect-other bi-allelic pairs (ambiguous)
    multiallelic_ambiguous_list = []  # Multiallelic sites that can't be uniquely resolved into a bi-allelic pair (multiallelic ambiguous)

    # Iterate over the chromosomes to select sites from and track failed edge cases
    for chrom in chromosomes_score:

        # Load the town's genomes for this chromosome
        filepath = f"data/genomes3/processed_chr{chrom}.csv"
        df = pd.read_csv(filepath)

        # Filter the polygenic score to this chromosome (t for trait)
        dt = dft.loc[dft.chrom == chrom].copy()

        # Merge the polygenic score with the town's genomes
        # Inner join ensures we only keep variants present in both datasets
        dt = pd.merge(dt, df, on="pos", how="inner")

        # Separate out multiallelic sites
        # These are identifed as having more than one ref, alt pair in the town's genomes
        multiallelic_sites = set(dt.groupby("pos").size()[lambda x: x > 1].index)

        # Put the multiallelic sites in a separate dataframe for processing
        dt_multiallelic = dt[dt.pos.isin(multiallelic_sites)]

        # Remove multiallelic sites from the polygenic chromosome dataframe
        dt = dt[~dt.pos.isin(multiallelic_sites)]

        # Get the sites where the effect allele is not the ref or alt allele (effect allele is missing)
        missing_sites = set(dt[(dt.effect_allele != dt.ref) & (dt.effect_allele != dt.alt)].pos)

        # Put the missing sites in a separate dataframe
        dt_missing = dt[dt.pos.isin(missing_sites)]

        # Remove the missing sites from the polygenic chromosome dataframe
        dt = dt[~dt.pos.isin(missing_sites)]

        # Get the sites where the effect allele is the ref allele and the other allele is not the alt (missing)
        effect_is_ref_other_is_missing_sites = set(dt[(dt.effect_allele == dt.ref) & (dt.other_allele != dt.alt)].pos)

        # Get the sites where the effect allele is the alt allele and the other allele is not the ref (missing)
        effect_is_alt_other_is_missing_sites = set(dt[(dt.effect_allele == dt.alt) & (dt.other_allele != dt.ref)].pos)

        # Ambiguous sites are where our ref-alt allele pair does not match the effect-other allele pair
        ambiguous_sites = set.union(effect_is_ref_other_is_missing_sites, effect_is_alt_other_is_missing_sites)

        # Put the ambiguous sites in a separate dataframe
        dt_ambiguous = dt[dt.pos.isin(ambiguous_sites)]

        # Remove the ambiguous sites from the polygenic chromosome dataframe
        dt = dt[~dt.pos.isin(ambiguous_sites)]

        # Iterate over the multiallelic sites to try resolving them
        multiallelic_sites = list(multiallelic_sites)

        # Create lists to store dataframes during iteration
        dt_resolved_list = []
        ds_multiallelic_ambiguous_list = []

        for site in multiallelic_sites:
            # Make a dataframe of the specific site 
            ds = dt_multiallelic[dt_multiallelic.pos == site].copy()

            # Select rows where the effect=other allele matches a ref-alt pair (could be either order)
            # This attempts to find a unique matching allele pair between the score and genome data
            ds_select = ds[((ds.effect_allele == ds.ref) & (ds.other_allele == ds.alt)) | ((ds.effect_allele == ds.alt) & (ds.other_allele == ds.ref))]

            # If the multiallelic site can be uniquely resolved to a pair, add it to the resolved list
            if len(ds_select) == 1:
                dt_resolved_list.append(ds_select)
            
            # If the multiallelic site can't be uniquely resolved, add it to the multiallelic ambiguous list
            else:
                ds_multiallelic_ambiguous_list.append(ds)

        # Concatenate the resolved sites to the main dataframe
        if len(dt_resolved_list) > 0:
            dt = pd.concat([dt] + dt_resolved_list)
        
        # Create the multiallelic ambiguous dataframe
        if len(ds_multiallelic_ambiguous_list) > 0:
            ds_multiallelic_ambiguous = pd.concat(ds_multiallelic_ambiguous_list)
        else:
            ds_multiallelic_ambiguous = pd.DataFrame()

        # Add a chrom column to all dataframes and put it as the first column
        dt["chrom"] = chrom
        dt = place_first(dt, "chrom")
        dt_missing["chrom"] = chrom
        dt_missing = place_first(dt_missing, "chrom")
        dt_ambiguous["chrom"] = chrom
        dt_ambiguous = place_first(dt_ambiguous, "chrom")
        ds_multiallelic_ambiguous["chrom"] = chrom
        ds_multiallelic_ambiguous = place_first(ds_multiallelic_ambiguous, "chrom")

        # Append non-empty dataframes to their respective lists
        if dt.shape[0] > 0:
            dts_list.append(dt)
        if dt_missing.shape[0] > 0:
            missing_list.append(dt_missing)
        if dt_ambiguous.shape[0] > 0:
            ambiguous_list.append(dt_ambiguous)
        if ds_multiallelic_ambiguous.shape[0] > 0:
            multiallelic_ambiguous_list.append(ds_multiallelic_ambiguous)

    # Concatenate all dataframes from their lists and reset indices
    if len(dts_list) > 0:
        dts = pd.concat(dts_list, ignore_index=True)
    else:
        dts = pd.DataFrame()
    if len(missing_list) > 0:
        dft_missing = pd.concat(missing_list, ignore_index=True)
    else:
        dft_missing = pd.DataFrame()
    if len(ambiguous_list) > 0:
        dft_ambiguous = pd.concat(ambiguous_list, ignore_index=True)
    else:
        dft_ambiguous = pd.DataFrame()
    if len(multiallelic_ambiguous_list) > 0:
        dft_multiallelic_ambiguous = pd.concat(multiallelic_ambiguous_list, ignore_index=True)
    else:
        dft_multiallelic_ambiguous = pd.DataFrame()

    # Save processed data to CSV files
    output_dir = f"data/score-genomes/{score_name}"
    ensure_dir(output_dir)

    if dts.shape[0] > 0:
        output_file = f"{output_dir}/selected.csv"
        dts.to_csv(output_file, index=False)

    # Save missing variants
    if dft_missing.shape[0] > 0:
        output_file = f"{output_dir}/missing.csv"
        dft_missing.to_csv(output_file, index=False)

    # Save ambiguous variants 
    if dft_ambiguous.shape[0] > 0:
        output_file = f"{output_dir}/ambiguous.csv"
        dft_ambiguous.to_csv(output_file, index=False)

    # Save multiallelic ambiguous variants
    if dft_multiallelic_ambiguous.shape[0] > 0:
        output_file = f"{output_dir}/multiallelic_ambiguous.csv" 
        dft_multiallelic_ambiguous.to_csv(output_file, index=False)

    # Fix strand flips and save
    df_strand_flipped = pd.DataFrame()

    if dft_missing.shape[0] > 0:
        df_strand_flipped = pd.concat([df_strand_flipped, dft_missing], ignore_index=True)
        
    if dft_multiallelic_ambiguous.shape[0] > 0:
        df_strand_flipped = pd.concat([df_strand_flipped, dft_multiallelic_ambiguous], ignore_index=True)

    if df_strand_flipped.shape[0] > 0:
        df_strand_flipped = fix_strand_flips(df_strand_flipped)
        output_file = f"{output_dir}/strand_flipped.csv"
        df_strand_flipped.to_csv(output_file, index=False)

    # Log match numbers
    log_match_numbers(score_name)


def fix_strand_flips(df):
    # Create new columns for analysis
    df = df.copy()
    
    def _get_complement(allele):
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        return complement.get(allele, allele)

    # Create complement alleles
    df['effect_complement'] = df['effect_allele'].apply(_get_complement)
    df['other_complement'] = df['other_allele'].apply(_get_complement)
    
    # Check if we need to flip (if complement matches better than original)
    needs_flip = (
        ((df['ref'] == df['effect_complement']) & (df['alt'] == df['other_complement'])) |
        ((df['ref'] == df['other_complement']) & (df['alt'] == df['effect_complement']))
    )
    
    # Where flip is needed, replace with complements
    df.loc[needs_flip, 'effect_allele'] = df.loc[needs_flip, 'effect_complement']
    df.loc[needs_flip, 'other_allele'] = df.loc[needs_flip, 'other_complement']
    
    # Drop temporary columns
    df = df.drop(['effect_complement', 'other_complement'], axis=1)
    
    return df[needs_flip]



def test_final_score_dataframe(dft: pd.DataFrame) -> None:
    """
    Implement basic QC checks of our score-town file
    """
    # Check that there is a bi-allelic match at every site between effect-other and alt-ref
    assert (((dft.effect_allele == dft.ref) & (dft.other_allele == dft.alt)) | ((dft.effect_allele == dft.alt) & (dft.other_allele == dft.ref))).all()

    # Assert that each site is unique
    assert one_per(dft, ["chrom", "pos"])
    assert one_per(dft, ["chrom", "pos", "effect_allele", "other_allele"])



def load_score_genome(score_name: str) -> dict:
    """
    Load a processed polygenic score from CSV files.
    Returns a dictionary containing available dataframes.
    """
    score_dir = f"data/score-genomes/{score_name}"
    files = {
        'selected': f"{score_dir}/selected.csv",
        'missing': f"{score_dir}/missing.csv", 
        'ambiguous': f"{score_dir}/ambiguous.csv",
        'multiallelic_ambiguous': f"{score_dir}/multiallelic_ambiguous.csv",
        'strand_flipped': f"{score_dir}/strand_flipped.csv",
        'original': f"data/scores/{score_name}.csv"
    }

    score_data = {}
    for name, filepath in files.items():
        if os.path.exists(filepath):
            score_data[name] = pd.read_csv(filepath)
        else:
            score_data[name] = pd.DataFrame()
            
    dft = score_data["selected"].copy()
    dft_strand_flipped = score_data["strand_flipped"].copy()
    
    # If there are strand flipped variants, add them to the final score dataframe
    if not dft_strand_flipped.empty:
        dft = pd.concat([dft, dft_strand_flipped], ignore_index=True)
    
    # Convert chrom to string so that it can be merged on
    dft.chrom = dft.chrom.astype(str)

    # Make the final score dataframe the selected + strand flipped
    score_data["final"] = dft

    # Implement basic QC checks of our score-town file
    test_final_score_dataframe(dft)

    return score_data

def log_match_numbers(score_name: str) -> None:
    """
    Compute and log the match rate between a polygenic score's SNPs and our genome data.

    """
    # Load the score data which contains SNP information
    score_data = load_score_genome(score_name)
    
    # Create dictionary to store chromosome-position sets for each SNP group
    chrom_pos_data = {}
    for group, dft in score_data.items():
        chrom_pos_data[group] = set(get_chrom_pos(dft)) if dft.shape[0] > 0 else set()
    
    # Combine all processed SNPs (selected, missing, ambiguous etc)
    chrom_pos_final = chrom_pos_data['final'] 
    original = chrom_pos_data['original']
    strand_flipped = chrom_pos_data['strand_flipped']

    # Create and save match statistics
    match_numbers = {"original": len(original), "final": len(chrom_pos_final), 
                     "strand_flipped": len(strand_flipped), 
                     "match_rate": len(chrom_pos_final) / len(original)}
    ensure_dir(f"data/logging/{score_name}")
    save_dict(match_numbers, f"data/logging/{score_name}/match_numbers.json")



def compute_score(score_name):
    """
    Compute polygenic scores for the european samples.
    
    Run time: about 1 second
    
    This function:
    1. Loads the processed score variants and their effect weights
    2. Validates the variant matching between score and genome data
    3. Computes weighted sums of alleles for each chromosome copy
    4. Combines maternal and paternal chromosome scores
    5. Maps scores to sample metadata and saves results
    
    The output dataframe contains:
    - name: Sample identifier
    - sex: Sex of the individual 
    - population_code: Population code (e.g. CEU)
    - child: Boolean indicating if sample is a child
    - father_name: Father's sample ID if known
    - mother_name: Mother's sample ID if known
    - population_name: Full population name
    - score: Computed polygenic score

    Args:
        score_name: Name of the polygenic score to compute

    Output is saved to data/computed/{score_name}/scores.csv
    """

    if os.path.exists(f"data/computed/{score_name}.csv"):
        return

    # Prepare for logging the score
    ensure_dir(f"data/logging/{score_name}")

    # Load processed score variants (combine the original selected and the strand-flipped)
    score_data = load_score_genome(score_name)

    # Get the final score dataframe
    dft = score_data["final"]

    # Check that the score does not contain X chromosome variants
    if "X" in set(dft.chrom):
        print(f"WARNING: Score {score_name} contains X chromosome variants, which are not supported yet. Skipping...")
        dft = dft[dft.chrom != "X"]
        dft.reset_index(drop=True, inplace=True)
        
    # The effect allele is either the ref or alt, demarcate the two cases
    effect_is_ref = dft.effect_allele == dft.ref
    effect_is_alt = dft.effect_allele == dft.alt

    # Count the number of sites where the effect allele is ref vs alt
    n_effect_is_ref = sum(effect_is_ref)
    n_effect_is_alt = sum(effect_is_alt)

    # Load genome data
    genome_names = load_genome_names()
    n_chromosomes = len(genome_names)
    european_sample_names = load_sample_names()
    n_samples = len(european_sample_names)
    assert n_chromosomes == n_samples * 2

    # Get (0, 1) boolean of whether the site has the alt allele for each chromosome in the genomes set (town) (haploid)
    dfn = dft[genome_names]

    # Get the weights of the effect allele (the weight per effect allele)
    weight = dft.effect_weight

    # Handle sites where effect allele is alt, if any exist
    if n_effect_is_alt > 0:
        # For the haploid samples
        # Get the number of effect alleles (0 or 1) - for sites where the effect allele is alt, this is just the number of alt alleles
        has_effect_allele_effect_is_alt = dfn[effect_is_alt].values
        assert set.issubset(set(has_effect_allele_effect_is_alt.flatten()), {0, 1})
        assert has_effect_allele_effect_is_alt.shape == (n_effect_is_alt, n_chromosomes)
        
        # Get the weights of the effect allele - for sites where the effect allele is alt
        weights_effect_is_alt = weight[effect_is_alt].values.reshape(-1,1)
        assert weights_effect_is_alt.shape == (n_effect_is_alt, 1)
        
        # Compute the loading for each site for each chromosome - for sites where the effect allele is alt
        loadings_effect_is_alt = weights_effect_is_alt * has_effect_allele_effect_is_alt
        assert loadings_effect_is_alt.shape == (n_effect_is_alt, n_chromosomes)

    # Demonstrate the flipping functions
    assert flip_zero_one(0) == 1
    assert flip_zero_one(1) == 0

    # Handle sites where effect allele is ref, if any exist
    if n_effect_is_ref > 0:
        # For the haploid samples
        # Get the number of effect alleles (0 or 1) - for sites where the effect allele is ref
        # We need to do flip ones to zeros and vice versa because 0 alt -> 1 ref -> 1 effect, and, 1 alt -> 0 ref -> 0 effect
        has_effect_allele_effect_is_ref = flip_zero_one(dfn[effect_is_ref].values)
        assert has_effect_allele_effect_is_ref.shape == (n_effect_is_ref, n_chromosomes)

        # Get the weights of the effect allele - for sites where the effect allele is ref
        weights_effect_is_ref = weight[effect_is_ref].values.reshape(-1,1)
        assert weights_effect_is_ref.shape == (n_effect_is_ref, 1)

        # Compute the loading for each site for each chromosome - for sites where the effect allele is ref
        loadings_effect_is_ref = weights_effect_is_ref * has_effect_allele_effect_is_ref
        assert loadings_effect_is_ref.shape == (n_effect_is_ref, n_chromosomes)

    # Sum up scores for each chromosome in the sample
    if n_effect_is_alt > 0 and n_effect_is_ref > 0:
        scores_chromosomes = loadings_effect_is_alt.sum(axis=0) + loadings_effect_is_ref.sum(axis=0)
    elif n_effect_is_alt > 0:
        scores_chromosomes = loadings_effect_is_alt.sum(axis=0)
    elif n_effect_is_ref > 0:
        scores_chromosomes = loadings_effect_is_ref.sum(axis=0)
    else:
        raise ValueError("No effect alleles found in the score data.")
    assert scores_chromosomes.shape == (n_chromosomes,)

    # Combine maternal and paternal chromosome scores (left chromosomes for each sample + right chromosomes for each sample)
    scores = scores_chromosomes[:n_samples] + scores_chromosomes[n_samples:]
    
    # Map scores to sample names 
    name_to_score = {}
    for name, score in zip(european_sample_names, scores):
        name_to_score[name] = score

    # Load and update sample metadata with scores
    dfs = load_df("european_samples")
    dfs["score"] = dfs.name.map(name_to_score)

    # Save processed data sample name and score to CSV files
    output_dir = f"data/computed"
    ensure_dir(output_dir)
    output_file = f"{output_dir}/{score_name}.csv"
    dfs.to_csv(output_file, index=False)


def characterize_score(score_name):
    """
    Characterize the score by computing the center and scale statistics using 5-fold CV.
    """

    if os.path.exists(f"data/characterized/{score_name}.csv"):
        return

    seed = 0

    # Define path to scores file
    filepath = f"data/computed/{score_name}.csv"

    # Read scores data into DataFrame
    dfs = pd.read_csv(filepath)

    # Filter for unrelated individuals
    dfu = dfs.loc[~dfs.related].copy()

    # Randomly shuffle the DataFrame and reset index
    dfu = dfu.sample(frac=1, random_state=seed).reset_index(drop=True)

    # Create 5-fold cross validation splits (0-4)
    n_folds = 5
    dfu['sample_fold'] = np.arange(len(dfu)) % n_folds

    # Iterate over the folds, robustly computing the center and scale for each
    entries = []
    for sample_fold in range(1, n_folds+1):
        dfut = dfu[dfu.sample_fold != sample_fold]
        
        median, sd = get_mad_normal(dfut.score)

        # sns.ecdfplot(data=dfut, x="score")
        # x = np.linspace(dfut.score.min(), dfut.score.max(), 1000)
        # y = norm.cdf(x, median, sd)
        # plt.plot(x, y, color="red", linewidth=2)

        # Perform Kolmogorov-Smirnov test on the normalized scores
        ks_pvalue = kstest((dfut.score - median) / sd, norm.cdf, alternative="two-sided").pvalue
        
        # Perform Shapiro-Wilk test on the normalized scores
        normalized_scores = (dfut.score - median) / sd
        _, shapiro_pvalue = shapiro(normalized_scores)

        entry = {"sample_fold": sample_fold, "center": median, "scale": sd, "ks_pvalue": ks_pvalue, "shapiro_pvalue": shapiro_pvalue}
        entries.append(entry)

    # Create a DataFrame for the current score
    df_score_stats_folds = pd.DataFrame(entries)

    # Save processed data sample name and score to CSV files
    output_dir = f"data/characterized"
    ensure_dir(output_dir)
    output_file = f"{output_dir}/{score_name}.csv"
    df_score_stats_folds.to_csv(output_file, index=False)


def obtain_score_transformation(score_name):
    """
    Obtains and saves the final score transformation file that will be applied to genomes.
    
    This function processes a score genome file and extracts the final transformation dataframe
    containing information like chromosome, alleles, effect weights, and variant positions.
    The output is saved as a CSV file that can be used to transform genome data.
    
    Example output columns:
    - chrom: Chromosome number
    - effect_allele/other_allele: The allele pairs for each variant
    - effect_weight: The weight to apply for the effect allele
    - rsid: Reference SNP ID
    - pos: Genomic position
    - cpra: Chromosome-Position-Reference-Alternate allele identifier
    - ref/alt: Reference and alternate alleles
    - acu/acr: Allele counts
    - cpra_match: Whether CPRA identifiers match
    
    Args:
        score_name: Name of the score to process
    """
    # Check if transformation file already exists
    if os.path.exists(f"data/scores2/{score_name}.csv"):
        return
    
    # Load the score genome data
    score_data = load_score_genome(score_name)

    # Load processed score variants (combine the original selected and the strand-flipped)
    score_data = load_score_genome(score_name)

    # Get the final score dataframe with transformation information
    dft = score_data["final"]

    # Added: 2025-04-10
    dft.pos = dft.pos.astype(int)

    # Filter out genome-specific columns to keep only transformation columns
    genome_names = load_genome_names()
    columns = [col for col in dft.columns if col not in genome_names]
    dft = dft[columns]
    
    # Save the transformation file
    output_dir = f"data/scores2"
    ensure_dir(output_dir)
    output_file = f"{output_dir}/{score_name}.csv"
    dft.to_csv(output_file, index=False)


def extract_european_genomes():
    chromosomes = [str(i) for i in range(1, 23)] + ["X"]
    for chrom in tqdm(chromosomes):
        extract_european_genome(chrom)
    os.system("rm -rf data/genomes/")


def extract_genotype_matrices():
    chromosomes = [str(i) for i in range(1, 23)]
    for chrom in tqdm(chromosomes):
        extract_genotype_matrix(chrom)
    extract_x_chr_genotype_matrix()
    os.system("rm -rf data/genomes2/")


def do_catalog_scores(catalog_score_names):
    
    print("Downloading catalog scores...")
    for score_name in tqdm(catalog_score_names):
        download_catalog_score(score_name)

    print("Processing catalog scores...")
    for score_name in tqdm(catalog_score_names):
        process_catalog_score(score_name)

    do_scores(catalog_score_names)


def do_scores(score_names):

    print("Intersecting scores...")
    for score_name in tqdm(score_names):
        try:
            intersect_score(score_name)
        except Exception as e:
            print(f"Error intersecting score {score_name}: {e}")

    print("Computing scores...")
    for score_name in tqdm(score_names):
        try:
            compute_score(score_name)
        except Exception as e:
            print(f"Error computing score {score_name}: {e}")

    print("Characterizing scores...")
    for score_name in tqdm(score_names):
        try:
            characterize_score(score_name)
        except Exception as e:
            print(f"Error characterizing score {score_name}: {e}")

    print("Obtaining score transformations...")
    for score_name in tqdm(score_names):
        try:
            obtain_score_transformation(score_name)
        except Exception as e:
            print(f"Error obtaining score transformation {score_name}: {e}")



def check_european_samples():
    """
    Performs validation checks on the European samples dataset.
    
    Loads the processed European samples dataframe and verifies:
    - The correct number of children (106)
    - Each child has both parents in the dataset
    - All parent IDs exist as samples in the dataset
    - The total number of people (532)
    - The correct number of unrelated individuals (214)
    
    Raises AssertionError if any checks fail.
    """
    dfs = load_df("european_samples")

    # Check number of children (should be 106)
    n_children = dfs.child.sum()
    print(f"Number of children: {n_children}")
    assert n_children == 106

    # Check that each child has both parents in dataset
    children = dfs[dfs.child == True]
    assert all(children.father_name.notna())
    assert all(children.mother_name.notna())

    # Check that all parents exist in dataset
    all_parents = set(children.father_name) | set(children.mother_name)
    all_names = set(dfs.name)
    assert all_parents.issubset(all_names)

    # Check total number of people (should be 532)
    total_people = len(dfs)
    print(f"Total people: {total_people}")
    assert total_people == 532

    # Check number of unrelated individuals (should be 212)
    n_in_trios = n_children * 3  # 106 children + 212 parents = 318
    n_unrelated = total_people - n_in_trios
    print(f"Number of unrelated individuals: {n_unrelated}")
    assert n_unrelated == 532 - (106 * 3) == 214

    print("All checks passed!")




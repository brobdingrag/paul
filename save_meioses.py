from inventory import *
from liftover import get_lifter
from consts import *

# Get the lifter
converter_forward = get_lifter('hg19', 'hg38', one_based=True)

# Open the pedigree simulator 
with open("ped-sim/father_child_bp.bp", "r") as f:
    bp_content = f.read()

# Split content into lines
lines = bp_content.strip().split('\n')

# Lists to store parsed data
records = []

print("Number of lines", len(lines))
for line in tqdm(lines):
    # Split line into basic components
    parts = line.split()
    if len(parts) < 4:  # Skip invalid lines
        continue
        
    label = parts[0]
    sex = parts[1]
    haplotype = parts[2]
    
    # Process chromosome segments
    current_pos = 3
    while current_pos < len(parts):
        segment = parts[current_pos]
        
        # Parse initial chromosome position
        if '|' in segment:
            chrom, start_pos = segment.split('|')
            start_pos = int(start_pos)
            current_founder = None
        # Parse breakpoint
        else:
            founder, end_pos = segment.split(':')
            end_pos = int(end_pos)
            
            # Create record for this segment
            records.append({
                'label': label,
                'sex': 'male' if sex == 's0' else 'female',
                'haplotype': haplotype,
                'chromosome': chrom,
                'start_position': start_pos,
                'end_position': end_pos,
                'founder_haplotype': founder
            })
            
            # Update start position for next segment
            start_pos = end_pos + 1
            
        current_pos += 1
        
# Create DataFrame
df = pd.DataFrame(records)

# Get the generation
df['generation'] = df['label'].str.extract(r'_g(\d+)-')[0].astype(int)

assert is_all(df[df.generation == 1].groupby("chromosome").end_position.nunique(), 1)
assert is_all(df[df.generation == 1].groupby("chromosome").start_position.nunique(), 1)

chrom_to_start = df[df.generation == 1].groupby("chromosome").start_position.first().to_dict()
chrom_to_end = df[df.generation == 1].groupby("chromosome").end_position.first().to_dict()

n_each_founder_haplotype = df[df.generation == 1].groupby(['label', 'haplotype']).size()

# Remove the first generation
df = df.loc[df.generation != 1]

# Ensure that the haplotypes are correct
assert set(df.groupby("haplotype").founder_haplotype.unique()["h1"]) == {"0", "1"}
assert set(df.groupby("haplotype").founder_haplotype.unique()["h0"]) == {"2", "3"}

# Get the individual number of the offspring
df['individual'] = df['label'].str.extract(r'i(\d+)')[0].astype(int)

n_breaks_dad = df[df.haplotype == "h0"].groupby("individual").size().quantile(0.5) - 22
n_breaks_mom = df[df.haplotype == "h1"].groupby("individual").size().quantile(0.5) - 22

# Dad's meiosis map for haplotype == h0 and Mom's meiosis map for haplotype == h1
df['parent'] = df.haplotype.map({"h0": "Dad", "h1": "Mom"})

# Get the grandparent
df['grandparent'] = df.founder_haplotype.map({"0": "Grandpa", "1": "Grandma", "2": "Grandpa", "3": "Grandma"})

# Drop haplotype columns as now redundant with parent and grandparent columns
df.drop(columns=['haplotype', 'founder_haplotype', 'label'], inplace=True)

# Rename some columns to be simpler and more descriptive
df.rename(columns={"chromosome": "chrom", "start_position": "pos_start_hg19", "end_position": "pos_end_hg19"}, inplace=True)

# Mark the start and end of the chromosome
df["start"] = False
df["end"] = False

# Iterate over the chromosomes
for chrom in set(df.chrom):
    pos_start = chrom_to_start[chrom]
    df.loc[(df.chrom == chrom) & (df.pos_start_hg19 == pos_start), "start"] = True
    df.loc[(df.chrom == chrom) & (df.pos_start_hg19 == pos_start), "pos_start_hg19"] = 0

    pos_end = chrom_to_end[chrom]

    df.loc[(df.chrom == chrom) & (df.pos_end_hg19 == pos_end), "end"] = True
    

def lift_forward_with_fallback(chrom: str, pos: int, window_size: int = 800, max_attempts: int = 600):
    """
    Attempts to lift coordinates forward, with fallback to nearby positions if initial position fails.
    
    Args:
        chrom: Chromosome number/name
        pos: Position to lift
        window_size: Size of window to search upstream/downstream
        max_attempts: Maximum number of positions to try before giving up
    
    Returns:
        Tuple of (original_position, lifted_position) or None if no valid lift found
    """
    chrom = str(chrom)
    
    def is_primary_assembly(result):
        """Check if result is on primary assembly and return position if it is."""
        if not result:
            return None
        for entry in result:
            if entry[0] == f"chr{chrom}":  # Only return primary assembly matches
                return entry[1]
        return None
    
    # Try original position first
    result = converter_forward.convert_coordinate(chrom, pos)
    primary_pos = is_primary_assembly(result)
    if primary_pos:
        return (pos, primary_pos)
    
    # Try positions up and downstream
    for i in range(1, max_attempts + 1):
        # Try upstream
        upstream_pos = pos - (i * window_size)
        result = converter_forward.convert_coordinate(chrom, upstream_pos)
        primary_pos = is_primary_assembly(result)
        if primary_pos:
            return (upstream_pos, primary_pos)
            
        # Try downstream
        downstream_pos = pos + (i * window_size)
        result = converter_forward.convert_coordinate(chrom, downstream_pos)
        primary_pos = is_primary_assembly(result)
        if primary_pos:
            return (downstream_pos, primary_pos)
            
    return None

# Start position
lifted_positions = df.apply(lambda row: lift_forward_with_fallback(row.chrom, row.pos_start_hg19) if not row.start else None, axis=1)
df['pos_start_hg19_query'] = lifted_positions.map(lambda x: x[0] if x else None)
df['pos_start_hg38'] = lifted_positions.map(lambda x: x[1] if x else None)

df.loc[df.start, "pos_start_hg19_query"] = 0
df.loc[df.start, "pos_start_hg38"] = 0

# End position
lifted_positions = df.apply(lambda row: lift_forward_with_fallback(row.chrom, row.pos_end_hg19) if not row.end else None, axis=1)
df['pos_end_hg19_query'] = lifted_positions.map(lambda x: x[0] if x else None)
df['pos_end_hg38'] = lifted_positions.map(lambda x: x[1] if x else None)

# Iterate over the chromosomes
for chrom in set(df.chrom):
    df.loc[(df.chrom == chrom) & (df.end), "pos_end_hg38"] = CHROMOSOME_LENGTHS[chrom]

# Failed individuals
failed_individuals = set(df[df.pos_end_hg38.isna()].individual).union(set(df[df.pos_start_hg38.isna()].individual))
n_failed = len(failed_individuals)

print("Number of failed individuals", n_failed)

df = df.loc[~df.individual.isin(failed_individuals)]

df.pos_start_hg38 = df.pos_start_hg38.astype(int)
df.pos_end_hg38 = df.pos_end_hg38.astype(int)

# Remove gaps in individual numbers
assert (df.individual.diff().iloc[1:] >= 0).all()
individuals_current = df.individual.drop_duplicates().to_numpy()
individuals_new = np.arange(n_unique(df.individual)) + 1
current_to_new = dict(zip(individuals_current, individuals_new))
df.rename(columns={"individual": "individual_2"}, inplace=True)
df["individual"] = df.individual_2.map(current_to_new)
assert n_unique(df.individual) == max(df.individual)

# Remove the previous individual column
df.drop(columns=["individual_2"], inplace=True)

# Remove the generation column 
assert is_all(df.generation, 2)
df.drop(columns=["generation"], inplace=True)

# Save the meioses
save_df(df, "meioses")


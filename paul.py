from consts import *
from inventory import *
from functools import cache
import multiprocessing as mp
from matplotlib.lines import Line2D
from scipy.stats import poisson
import urllib.request
from liftover import get_lifter


# Get the lifter
converter_forward = get_lifter('hg19', 'hg38', one_based=True)


def iterate_autosomes():
    for chrom in range(1, 23):
        yield chrom


def get_couples_embryos(generation, couple_num):

    dc = load_couples(generation)

    couple = dc[dc.couple_num == couple_num].squeeze()

    check_parents_sex(generation, couple.name_dad, couple.name_mom)

    if couple.n_embryos == 0:
        return

    ensure_dir(f"data/embryo_score_genomes/{generation}/{couple_num}")
    meiosis_nums = {}
    for embryo_num in range(1, couple.n_embryos + 1):

        sperm, sperm_meiosis_num = get_gamete_states(generation, couple.name_dad)
        egg, egg_meiosis_num = get_gamete_states(generation, couple.name_mom)

        meiosis_nums[embryo_num] = {
            "sperm_meiosis_num": sperm_meiosis_num,
            "egg_meiosis_num": egg_meiosis_num,
        }

        dft = sperm.merge(egg, on=["chrom", "pos"]) 
        save_df(dft, f"embryo_score_genomes/{generation}/{couple_num}/{embryo_num}")

    ensure_dir(f"data/embryo_meiosis_nums/{generation}")
    save_dict(meiosis_nums, f"data/embryo_meiosis_nums/{generation}/{couple_num}")


def get_generation_embryos(generation, ncores=8):
    dc = load_couples(generation)
    
    # Prepare arguments for parallel processing
    tasks = [(generation, couple_num) for couple_num in dc.couple_num]

    with mp.Pool(processes=ncores) as pool:
        pool.starmap(get_couples_embryos, tasks)


def load_couples(generation):
    dc = load_df(f"data/couples/{generation}/couples")

    # Check that each dad and mom is only in one couple
    # (Check for monogamy)
    assert dc.name_dad.is_unique
    assert dc.name_mom.is_unique

    return dc


def get_couples(generation):

    # Get the sex of each person
    name_to_sex = load_name_to_sex(generation)

    # Transform sex to a dataframe
    dc = pd.Series(name_to_sex).reset_index().rename(columns={"index": "name", 0: "sex"})

    # Get the number of males and females
    n_females = dc[dc.sex == "female"].shape[0]
    n_males = dc[dc.sex == "male"].shape[0]

    # The number of couples is the minimum of the number of males and females
    n_couples = min(n_females, n_males)

    # Sample the males and females
    potential_dads = dc[dc.sex == "male"].sample(n=n_couples).name.tolist()
    potential_moms = dc[dc.sex == "female"].sample(n=n_couples).name.tolist()

    # Keep track of the moms that have been paired
    paired_moms = set()

    # Initialize list of couples
    couples = []

    # Initialize the couple number counter
    couple_num = 1

    # Iterate over the potential dads
    for name_dad in potential_dads:

        # Iterate over the potential moms
        for name_mom in potential_moms:

            # Skip if the mom has already been paired
            if name_mom in paired_moms:
                continue

            # Skip if the dad and mom are related
            if are_related(name_dad, name_mom, generation):
                continue

            # Create a new couple
            couple = {
                "couple_num": couple_num,
                "name_dad": name_dad,
                "name_mom": name_mom,
                "n_embryos": get_n_embryos(generation)
            }

            # Add the couple to the list
            couples.append(couple)

            # Increment the couple number
            couple_num += 1

            # Add the mom to the paired moms set
            paired_moms.add(name_mom)

            break

    dc = pd.DataFrame(couples)

    # Save the couples
    ensure_dir(f"data/couples/{generation}")
    save_df(dc, f"data/couples/{generation}/couples")


def get_gamete_states(generation, name):
    """Get the states of the score sites in a person's (random) gamete"""
    # Get a random meiosis for the person
    dm, meiosis_num = get_meiosis(generation, name)

    # Load the states of the score sites in the person
    dft = load_score_genome(generation, name)

    # Determine the states of the gamete's score alleles
    # Initialize the parental chromosome of the gamete
    dft['parental_chrom'] = ""
    for chrom in iterate_autosomes():  
        dmt = dm[dm.chrom == chrom]
        for _, chunk in dmt.iterrows():
            parental_chrom = "paternal" if chunk.grandparent == "Grandpa" else "maternal"
            chunk_mask = (dft.chrom == chrom) & (dft.pos.between(chunk.pos_start, chunk.pos_end))
            dft.loc[chunk_mask, 'parental_chrom'] = parental_chrom

    assert set(dft.parental_chrom) == {"paternal", "maternal"}

    dft['gamete'] = np.nan
    dft.loc[dft.parental_chrom == "paternal", "gamete"] = dft['paternal']
    dft.loc[dft.parental_chrom == "maternal", "gamete"] = dft['maternal']
    dft.gamete = dft.gamete.astype(int)
    assert dft.gamete.isna().sum() == 0

    dft = dft[['chrom', 'pos', 'gamete']]
    
    sex = load_sex(generation, name)
    haploid_name = "paternal" if sex == "male" else "maternal"
    dft.rename(columns={"gamete": haploid_name}, inplace=True)

    return dft, meiosis_num


def load_pmf_n_embryos():
    lamb = N_EMBRYOS_PER_IVF * N_IVF_PER_COUPLE
    dp = pd.Series(poisson.pmf(k=range(33), mu=lamb))
    dp *= (1 - FRACTION_INFERTILE)
    dp.loc[0] += FRACTION_INFERTILE

    dp = dp.reset_index()
    dp = dp.rename(columns={'index': 'n_embryos', 0: 'proba'})
    return dp


def get_n_embryos(generation):
    """Randomly draw a number of viable embryos for a couple."""
    dp = load_pmf_n_embryos()
    n_embryos = np.random.choice(dp.n_embryos, p=dp.proba)
    while n_embryos == 0 and generation == 0:
        n_embryos = np.random.choice(dp.n_embryos, p=dp.proba)
    return n_embryos


def load_name_to_sex(generation):
    name_to_sex = load_dict(f"sexes/{generation}/name_to_sex")
    name_to_sex = {int(k): v for k, v in name_to_sex.items()}
    return name_to_sex


def load_sex(generation, name):
    name_to_sex = load_name_to_sex(generation)
    sex = name_to_sex[name]
    assert sex in ["male", "female"]
    return sex


def load_n_meioses(force_reload=False):
    if not os.path.exists("data/n_meioses.txt") or force_reload:
        dm = load_meioses()
        n_meioses = int(dm.meiosis_num.nunique())
        save_int(n_meioses, "n_meioses")
    n_meioses = load_int("n_meioses")
    return n_meioses


def load_meioses(force_reload=False):
    if not os.path.exists("data/meioses_fixed.csv") or force_reload:
        dm = load_df("meioses")
        dm = dm[dm.chrom != "X"]
        dm.chrom = dm.chrom.astype(int)
        exclude_cols = [
            "sex",
            'pos_start_hg19', 'pos_end_hg19', 
            'pos_start_hg19_query', 'pos_end_hg19_query'
            ]
        dm.rename(columns={
            "pos_start_hg38": "pos_start",
            "pos_end_hg38": "pos_end"
            }, inplace=True)
        dm = dm.drop(columns=exclude_cols)
        dm.rename(columns={"individual": "meiosis_num"}, inplace=True)
        dm.reset_index(drop=True, inplace=True)
        
        # Make the pos_start and pos_end next to each other
        # a dozen or so hg19 -> hg38 alignment issues
        for i, row in dm.iterrows():
            if i == len(dm) - 2:
                break
            row_next = dm.iloc[i+ 1]
            assert row.name == i
            assert row_next.name == i + 1
            if row.chrom == row_next.chrom and row.meiosis_num == row_next.meiosis_num:
                if row.pos_end + 1 != row_next.pos_start:
                    dm.loc[dm.index == i + 1, "pos_start"] = row.pos_end + 1
                    # difference = row_next.pos_start - (row.pos_end + 1)
                    # counts[difference] += 1
        
        # Check that all the chunks are adjacent now
        for i, row in dm.iterrows():
            if i == len(dm) - 2:
                break
            row_next = dm.iloc[i+ 1]
            assert row.name == i
            assert row_next.name == i + 1
            if row.chrom == row_next.chrom and row.meiosis_num == row_next.meiosis_num:
                assert row.pos_end + 1 == row_next.pos_start
        
        save_df(dm, "meioses_fixed")

    return load_df("meioses_fixed")





def get_random_meiosis():
    n_meioses = load_n_meioses()
    random_meiosis = np.random.choice(n_meioses) + 1
    return random_meiosis


def get_random_sex():
    sex = np.random.choice(["male", "female"])
    sex = str(sex)
    return sex


def get_meiosis(generation, name):
    dm = load_meioses()
    sex = load_sex(generation, name)
    parent = "Mom" if sex == "female" else "Dad"
    dm = dm[dm.parent == parent]

    meiosis_num = get_random_meiosis()
    dm = dm[dm.meiosis_num == meiosis_num]

    dm = dm.reset_index(drop=True)
    return dm, meiosis_num


def check_parents_sex(generation, name_dad, name_mom):
    sex_dad = load_sex(generation, name_dad)
    assert sex_dad == "male"
    sex_mom = load_sex(generation, name_mom)
    assert sex_mom == "female"


def load_score_genome(generation, name):
    dft = load_df(f"score_genomes/{generation}/{name}")
    return dft


def load_score_embryo(generation, couple_num, embryo_num):
    dft = load_df(f"embryo_score_genomes/{generation}/{couple_num}/{embryo_num}")
    return dft


def score_person(generation, name):
    dft = load_score_genome(generation, name)
    score = score_genome(dft)
    return score


def score_embryo(generation, couple_num, embryo_num):
    dft = load_score_embryo(generation, couple_num, embryo_num)
    embryo_score = score_genome(dft)
    return embryo_score


def score_genome(dft):
    """Score the genome of an individual (embryo or person)
    with these columns: ['chrom', 'pos', 'paternal', 'maternal']
    """
    dfs = load_df("height_score")
    dft = dft.merge(dfs, on=['chrom', 'pos'])

    # The effect allele is either the ref or alt, demarcate the two cases
    effect_is_ref = dft.effect_allele == dft.ref
    effect_is_alt = dft.effect_allele == dft.alt

    assert all(effect_is_ref | effect_is_alt)

    # Count the number of sites where the effect allele is ref vs alt
    n_effect_is_ref = sum(effect_is_ref)
    n_effect_is_alt = sum(effect_is_alt)

    # Load genome data
    n_chromosomes = 2
    chromosome_names = ['paternal', 'maternal']

    # Get (0, 1) boolean of whether the site has the alt allele for each chromosome
    dfn = dft[chromosome_names]

    # Get the weights of the effect allele (the weight per effect allele)
    weight = dft.effect_weight

    # Handle sites where effect allele is alt, if any exist
    if n_effect_is_alt > 0:

        # For the haploid samples
        # Get the number of effect alleles (0 or 1) - for sites where the effect allele is alt, 
        # this is just the number of alt alleles
        has_effect_allele_effect_is_alt = dfn[effect_is_alt].values
        assert set.issubset(set(has_effect_allele_effect_is_alt.flatten()), {0, 1})
        assert has_effect_allele_effect_is_alt.shape == (n_effect_is_alt, n_chromosomes)
        
        # Get the weights of the effect allele - for sites where the effect allele is alt
        weights_effect_is_alt = weight[effect_is_alt].values.reshape(-1,1)
        assert weights_effect_is_alt.shape == (n_effect_is_alt, 1)
        
        # Compute the loading for each site for each chromosome - for sites where the effect allele is alt
        loadings_effect_is_alt = weights_effect_is_alt * has_effect_allele_effect_is_alt
        assert loadings_effect_is_alt.shape == (n_effect_is_alt, n_chromosomes)

    def flip_zero_one(x):
        """Flips 0 to 1 and 1 to 0."""
        return (x - 1) * -1

    # Demonstrate the flipping functions
    assert flip_zero_one(0) == 1
    assert flip_zero_one(1) == 0

    # Handle sites where effect allele is ref, if any exist
    if n_effect_is_ref > 0:
        # For the haploid samples
        # Get the number of effect alleles (0 or 1) - for sites where the effect allele is ref
        # We need to do flip ones to zeros and vice versa because 0 alt -> 1 ref -> 1 effect, 
        # and, 1 alt -> 0 ref -> 0 effect
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

    score = sum(scores_chromosomes)

    return score


def get_couples_living_embryos(generation, couple_num):
    dc = load_couples(generation)
    couple = dc[dc.couple_num == couple_num].squeeze()

    if couple.n_embryos == 0:
        return []

    embryo_scores = {embryo_num: score_embryo(generation, couple_num, embryo_num) 
                    for embryo_num in range(1, couple.n_embryos + 1)}

    de = pd.Series(embryo_scores).to_frame().reset_index()
    de.rename(columns={"index": "embryo_num", 0: "score"}, inplace=True)
    de.sort_values("score", ascending=False, inplace=True)
    de['live_birth'] = [np.random.binomial(1, BIRTH_RATE) for _ in range(len(de))]
    dl = de[de.live_birth == 1].head(N_DESIRED_BIRTHS)
    living_embryos = dl.embryo_num.tolist()
    return living_embryos


def select_generation_embryos(generation):
    """
    Assumes I already created the couples and their embryos
    with get_couples(generation) and get_generation_embryos(generation)

    This does selection, birth, and formation of next generation.
    """

    dc = load_couples(generation)

    # Initialize the name counter
    name = 1

    # Initialize name to sex for next generation
    name_to_sex = {}

    # Iterate over the couples
    for couple_num in dc.couple_num:
        couple = dc[dc.couple_num == couple_num].squeeze()

        # Get the living embryos of the couple (selection)
        living_embryos = get_couples_living_embryos(generation, couple_num)
        
        # If the couple has no living embryos, skip
        if len(living_embryos) == 0:
            continue

        # Iterate over the living embryos
        for embryo_num in living_embryos:
            # Randomly assign sex
            sex = get_random_sex()

            # Store the parental data
            data = {
                "name": name,
                "sex": sex,
                "embryo_num": embryo_num,
                "dad": int(couple.name_dad),
                "mom": int(couple.name_mom),
            }

            # Store the sex
            name_to_sex[name] = sex

            # Save the parental data
            ensure_dir(f"data/parents/{generation + 1}/")
            save_dict(data, f"data/parents/{generation + 1}/{name}")

            # Save the embryo score genome as next generation person
            ensure_dir(f"data/score_genomes/{generation + 1}")
            embryo_path = f"data/embryo_score_genomes/{generation}/{couple_num}/{embryo_num}.csv"
            person_path = f"data/score_genomes/{generation + 1}/{name}.csv"
            command = f"cp {embryo_path} {person_path}"
            os.system(command)

            name += 1

    # Save the sex data for next generation
    ensure_dir(f"data/sexes/{generation + 1}/")
    save_dict(name_to_sex, f"data/sexes/{generation + 1}/name_to_sex")


def load_generation_names(generation):
    name_to_sex = load_name_to_sex(generation)
    return list(name_to_sex.keys())


def load_population_sizes():
    sizes = np.array([len(load_generation_names(i)) for i in range(N_GENERATIONS + 1)])
    return sizes


def are_related(name_1, name_2, generation, degree=1):
    if generation == 0:
        return False
    
    if degree >= DEGREE_UNRELATED:
        return False
    
    parents_1 = load_dict(f"data/parents/{generation}/{name_1}")
    parents_2 = load_dict(f"data/parents/{generation}/{name_2}")

    dad_1, mom_1 = parents_1['dad'], parents_1['mom']
    dad_2, mom_2 = parents_2['dad'], parents_2['mom']
    
    if dad_1 == dad_2 or mom_1 == mom_2:
        print('hi')
        return True
    
    dad_1_dad_2_related = are_related(dad_1, dad_2, generation - 1, degree + 1)
    
    if dad_1_dad_2_related:
        return True
    dad_1_mom_2_related = are_related(dad_1, mom_2, generation - 1, degree + 1)
    if dad_1_mom_2_related:
        return True
    mom_1_dad_2_related = are_related(mom_1, dad_2, generation - 1, degree + 1)
    if mom_1_dad_2_related:
        return True
    mom_1_mom_2_related = are_related(mom_1, mom_2, generation - 1, degree + 1)
    if mom_1_mom_2_related:
        return True

    return False


def get_intergenerational_scores(generation=0):
    """Gets the scores of the parents of the `generation` and their embryos."""

    dc = load_df(f"couples/{generation}/couples")
    dc.rename(columns={"name_dad": "dad", "name_mom": "mom"}, inplace=True)

    assert dc.mom.is_unique, "not monogamous"
    assert dc.dad.is_unique, "not monogamous"

    dc['dad_score'] = dc.dad.apply(lambda x: score_person(generation, x))
    dc['mom_score'] = dc.mom.apply(lambda x: score_person(generation, x))

    parent_data = []
    for kid_name in load_generation_names(generation + 1):
        parent_data.append(load_dict(f"parents/{generation + 1}/{kid_name}"))
    dp = pd.DataFrame(parent_data)

    data = []
    for _, couple in dc.iterrows():
        if couple.n_embryos == 0:
            continue
        meiosis_nums = load_dict(f"embryo_meiosis_nums/{generation}/{int(couple.couple_num)}")
        meiosis_nums = {int(k): v for k, v in meiosis_nums.items()}
        for embryo_num in range(1, int(couple.n_embryos) + 1):
            entry = couple.to_dict()
            entry['embryo_num'] = embryo_num
            entry['sperm_meiosis_num'] = meiosis_nums[embryo_num]['sperm_meiosis_num']
            entry['egg_meiosis_num'] = meiosis_nums[embryo_num]['egg_meiosis_num']
            data.append(entry)

    df = pd.DataFrame(data)

    df = df.astype({'couple_num': int, 'dad': int, 'mom': int, 'n_embryos': int})

    dp = dp.drop("sex", axis=1)
    embryo_ids = ['dad', 'mom', 'embryo_num']
    df = df.set_index(embryo_ids).join(dp.set_index(embryo_ids)).reset_index()
    df['selected'] = df.name.notna()

    df.name = df.name.fillna(-1).astype(int)
    df['embryo_score'] = df.apply(lambda row: score_embryo(
                        generation, row.couple_num, row.embryo_num), axis=1)
    df['dad_mom_avg_score'] = (df.dad_score + df.mom_score) / 2
    df['embryo_score_normalized'] = df.embryo_score - df.dad_mom_avg_score

    df['generation'] = generation

    df = place_first(df, 'generation')

    return df


def get_generation_scores(generation=0):
    names = load_generation_names(generation)
    name_to_score = {name: score_person(generation, name) for name in names}
    name_to_sex = load_name_to_sex(generation)
    df = pd.DataFrame(names, columns=['name'])
    df['score'] = df.name.map(name_to_score)
    df['sex'] = df.name.map(name_to_sex)
    df['generation'] = generation

    df = place_first(df, 'generation')
    return df


def save_intergenerational_scores():
    generations = list(range(N_GENERATIONS))
    
    with mp.Pool(processes=N_GENERATIONS) as pool:
        data = pool.map(get_intergenerational_scores, generations)
    
    df = pd.concat(data)
    save_df(df, "intergenerational_scores")


def save_generation_scores():
    generations = list(range(N_GENERATIONS + 1))
    
    with mp.Pool(processes=N_GENERATIONS) as pool:
        data = pool.map(get_generation_scores, generations)
    
    df = pd.concat(data)
    save_df(df, "generation_scores")


def load_intergenerational_scores():
    di = load_df("intergenerational_scores")
    return di


def load_generation_scores(generation=None):
    dg = load_df("generation_scores")
    if generation is not None:
        dg = dg[dg.generation == generation]
    return dg


def save_freq_vs_weight():
    dfs = load_df("height_score")
    dfs.set_index(['chrom', 'pos'], inplace=True)

    for generation in tqdm(range(N_GENERATIONS + 1)):
        
        names = load_generation_names(generation)
        
        dc = pd.Series(index=dfs.index, data=np.zeros(len(dfs)), dtype=int)
        
        for name in names:
            dft = load_df(f"score_genomes/{generation}/{name}")
            person_count = dft.set_index(['chrom', 'pos']).sum(axis=1)
            dc += person_count
        
        freq = dc / len(names) / 2
        
        dfs[f"freq_gen_{generation}"] = freq

    effect_is_ref = dfs.effect_allele == dfs.ref
    effect_is_alt = dfs.effect_allele == dfs.alt
    dft['alt_weight'] = np.nan
    dfs.loc[effect_is_alt, "alt_weight"] = dfs.effect_weight
    dfs.loc[effect_is_ref, "alt_weight"] = -1 * dfs.effect_weight
    dfs['minor_weight'] = np.nan
    alt_is_major_gen_0 = dfs.freq_gen_0 >= 0.5
    alt_is_minor_gen_0 = dfs.freq_gen_0 < 0.5
    dfs.loc[alt_is_minor_gen_0, "minor_weight"] = dfs.alt_weight
    dfs.loc[alt_is_major_gen_0, "minor_weight"] = -1 * dfs.alt_weight
    freq_cols = [f"freq_gen_{i}" for i in range(N_GENERATIONS + 1)]
    dfs.loc[alt_is_major_gen_0, freq_cols] = 1 - dfs.loc[alt_is_major_gen_0, freq_cols]

    save_df(dfs.reset_index(), "freq_vs_weight")



def load_freq_vs_weight():
    dfs = load_df("freq_vs_weight")
    dfs = dfs.set_index(['chrom', 'pos'])
    return dfs



def get_gamete_ancestry_chunks(parent_ancestry, parent_meiosis):
        
    gamete_ancestry_chunks = []

    for chrom in iterate_autosomes():
        ancestry_chunks = parent_ancestry[parent_ancestry.chrom == chrom]
        gamete_chunks = parent_meiosis[parent_meiosis.chrom == chrom]

        for _, gamete_chunk in gamete_chunks.iterrows():
            parental_chrom = "paternal" if gamete_chunk.grandparent == "Grandpa" else "maternal"
            overlapping_ancestry_chunks = ancestry_chunks[
                                            (ancestry_chunks.member == parental_chrom) & 
                                            (ancestry_chunks.pos_start <= gamete_chunk.pos_end) & 
                                            (ancestry_chunks.pos_end >= gamete_chunk.pos_start)
                                            ]
            for _, ancestry_chunk in overlapping_ancestry_chunks.iterrows():
                gamete_ancestry_chunk = {
                    "chrom": chrom,
                    "pos_start": max(ancestry_chunk.pos_start, gamete_chunk.pos_start),
                    "pos_end": min(ancestry_chunk.pos_end, gamete_chunk.pos_end),
                    "ancestor": ancestry_chunk.ancestor,
                    "ancestor_member": ancestry_chunk.ancestor_member,
                }
                gamete_ancestry_chunks.append(gamete_ancestry_chunk)
            
    return pd.DataFrame(gamete_ancestry_chunks)


def save_ancestry_chunks():
        
    di = load_intergenerational_scores() 
    dm = load_meioses()

    # Initialize the ancestry datafram
    generation = 0
    data = []
    for name in load_generation_names(generation):
        
        for chrom in iterate_autosomes():
            
            ancestor = name

            for member in ["paternal", "maternal"]:

                entry = {
                    "generation": generation, 
                    "name": name,
                    "chrom": chrom,
                    "pos_start": 0,
                    "pos_end": CHROMOSOME_LENGTHS[str(chrom)],
                    "member": member,
                    "ancestor": ancestor,
                    "ancestor_member": member,
                }

                data.append(entry)

    dfa = pd.DataFrame(data)

    for generation in range(1, N_GENERATIONS + 1):

        print(f"Generation {generation}")

        dc = load_couples(generation-1)

        names = load_generation_names(generation)

        data_generation = []
        for name in tqdm(names):
                
            # Get the genesis of the person 
            genesis = di[(di.generation == generation-1) & (di.name == name)].squeeze()

            # Get the dad and mom of the person
            dad, mom = genesis.dad, genesis.mom

            # Get the meiosis numbers of the dad and mom
            sperm_meiosis_num, egg_meiosis_num = genesis.sperm_meiosis_num, genesis.egg_meiosis_num

            # Double check that the dad and mom were a couple in the previous generation
            couple = dc[(dc.name_dad == dad) & (dc.name_mom == mom)].squeeze()
            assert not couple.empty

            # Get Dad's ancestry
            dad_ancestry = dfa[(dfa.generation == generation-1) & (dfa.name == dad)]

            # Get Dad's meiosis (his sperm)
            sperm_meiosis = dm[(dm.meiosis_num == sperm_meiosis_num) & (dm.parent == "Dad")]

            # Get Mom's ancestry
            mom_ancestry = dfa[(dfa.generation == generation-1) & (dfa.name == mom)]

            # Get Mom's meiosis (her egg)
            egg_meiosis = dm[(dm.meiosis_num == egg_meiosis_num) & (dm.parent == "Mom")]

            # Get the gamete ancestry chunks
            sperm_ancestry_chunks = get_gamete_ancestry_chunks(dad_ancestry, sperm_meiosis)
            egg_ancestry_chunks = get_gamete_ancestry_chunks(mom_ancestry, egg_meiosis)
            
            sperm_ancestry_chunks['member'] = "paternal"
            egg_ancestry_chunks['member'] = "maternal"

            dfa_person = pd.concat([sperm_ancestry_chunks, egg_ancestry_chunks])
            dfa_person['name'] = name

            data_generation.append(dfa_person)

        dfa_generation = pd.concat(data_generation)
        dfa_generation['generation'] = generation
        dfa = pd.concat([dfa, dfa_generation])

    save_df(dfa, "ancestry")

@cache
def load_ancestry(generation=None):
    dfa = load_df("ancestry")
    if generation is not None:
        dfa = dfa[dfa.generation == generation]
    return dfa



@cache
def load_name_to_id_1kg():
    id_1kg_to_name = load_id_1kg_to_name()
    name_to_id_1kg = {name: id_1kg for id_1kg, name in id_1kg_to_name.items()}
    return name_to_id_1kg


@cache
def load_id_1kg_to_name():
    dfs = load_df("european_samples")
    dfs = dfs[~dfs.related]
    id_1kg_to_name = {id_1kg: i+1 for i, id_1kg in enumerate(dfs.name)}
    return id_1kg_to_name





def store_initial_generation():
    generation = 0
    
    dfs = load_df("european_samples")
    dfs = dfs[~dfs.related]
    # rename_map = {name: i+1 for i, name in enumerate(dfs.name)}
    rename_map = load_id_1kg_to_name()

    df = load_df("height_states")

    ensure_dir(f"data/score_genomes/{generation}")
    for name_orig, name in rename_map.items():
        dt = df[['chrom', 'pos', f"{name_orig}_1", f"{name_orig}_2"]].copy()
        dt.rename(columns={f"{name_orig}_1": "paternal", f"{name_orig}_2": "maternal"}, inplace=True)
        save_df(dt, f"score_genomes/{generation}/{name}")
        
    name_orig_to_sex = dfs.set_index("name").sex.to_dict()
    name_to_sex = {rename_map[name]: name_orig_to_sex[name] for name in rename_map}
    ensure_dir(f"data/sexes/{generation}")
    save_dict(name_to_sex, f"sexes/{generation}/name_to_sex")





def store_height_score():
    df = load_df("height_states")
    cols = ['chrom', 'pos', 'effect_allele', 'ref', 'alt', 'other_allele', 'effect_weight', 'rsid']
    dfs = df[cols]
    save_df(dfs, "height_score")




def lift_row_forward(row, return_type='pos'):

    assert return_type in ['pos', 'strand', 'both']
    conversion = converter_forward.convert_coordinate(str(row.chrom), row.pos_hg19)
    if len(conversion) > 1:
        raise ValueError(f"Multiple conversions for {row.chrom}:{row.pos_hg19}")
    if len(conversion) == 0:
        if return_type == 'pos':
            return None
        elif return_type == 'strand':
            return None
        elif return_type == 'both':
            return None, None
    conversion = conversion[0]
    if f"chr{row.chrom}" != conversion[0]:
        pos, strand = None, None
    else:
        pos = conversion[1]
        strand = conversion[2]
    if return_type == 'pos':
        return pos
    elif return_type == 'strand':
        return strand
    elif return_type == 'both':
        return pos, strand


def store_ukb_pca_map():

    print("Downloading SNP PCA map...")

    file_url = "biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/snp_pca_map.txt"

    urllib.request.urlretrieve(f"https://{file_url}", "data/snp_pca_map.txt")

    # Load the SNP PCA mapping data from the UK Biobank
    dt = pd.read_csv(f"data/snp_pca_map.txt", sep="\t", header=None)

    # Rename columns based on the file structure
    # First 7 columns contain SNP information, followed by PC loadings
    column_names = ['chrom', 'rsid', 'cm_location', 'pos_hg19', 
                    'allele1', 'allele2', 'freq_allele1']

    dt = dt.iloc[:, :7]
    dt.columns = column_names

    # Remove cm_location column if it contains only one unique value (uninformative)
    if n_unique(dt.cm_location) == 1:
        dt.drop("cm_location", axis=1, inplace=True)

    print("Lifting SNPs to hg38...")
    # Convert hg19 positions to hg38 using the lift_row_forward function
    # This returns both the new position and strand information
    dt['pos'] = dt.apply(lambda row: lift_row_forward(row, return_type='pos'), axis=1)
    dt['strand'] = dt.apply(lambda row: lift_row_forward(row, return_type='strand'), axis=1)

    # Reorganize columns for better readability
    dt = place_first(dt, ['chrom', 'pos', 'strand', 'allele1', 'allele2', 'freq_allele1'])

    # Save SNPs that couldn't be lifted to hg38 (NA positions) for later analysis
    dtu = dt[dt.pos.isna()].copy()
    print(f"{len(dtu)} SNPs couldn't be lifted to hg38")

    # Continue with SNPs that were successfully lifted to hg38
    dt = dt[dt.pos.notna()]
    assert no_na(dt.pos)  # Verify no missing positions
    assert no_na(dt.strand)  # Verify no missing strand information
    dt.pos = dt.pos.astype(int)  # Convert positions to integers

    # Fix strand flips for SNPs on the negative strand
    # Convert alleles to their complementary bases
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    dt.loc[dt.strand == "-", "allele1"] = dt.loc[dt.strand == "-"]['allele1'].map(complement)
    dt.loc[dt.strand == "-", "allele2"] = dt.loc[dt.strand == "-"]['allele2'].map(complement)

    # Final formatting improvements
    dt.rename(columns={"allele1": "allele_1", "allele2": "allele_2",
                       "freq_allele1": "freq_allele_1"}, inplace=True)
    dt = place_first(dt, ['chrom', 'pos', 'strand'])
    
    save_df(dt, "ukb_pca_map")
    print("Done!")


def load_ukb_pca_map():
    dt = load_df("ukb_pca_map")
    return dt


def get_1kg_snps(chrom, positions_needed, chunksize=100_000):
    # Process the file in chunks
    chunks = []
    for chunk in pd.read_csv(f"data/genomes3/processed_chr{chrom}.csv", chunksize=chunksize, dtype={'pos': int}):
        filtered_chunk = chunk[chunk['pos'].isin(positions_needed)]
        chunks.append(filtered_chunk)

    # Combine into a single DataFrame
    dv = pd.concat(chunks)

    dv['chrom'] = chrom
    dv = place_first(dv, "chrom")

    return dv

def get_positions_needed(df, chrom):
    return set(df[df.chrom == chrom].pos)


def store_1kg_array(dt, array_name, ncores=8):
    tasks = [(chrom, get_positions_needed(dt, chrom)) for chrom in iterate_autosomes()]
        
    with mp.Pool(processes=ncores) as pool:
        results = pool.starmap(get_1kg_snps, tasks)

    dvt = pd.concat(results, ignore_index=True)

    dvt = dt.merge(dvt, on=['chrom', 'pos'])

    dvt['allele_1_complement'] = dvt.allele_1.map(get_complement)
    dvt['allele_2_complement'] = dvt.allele_2.map(get_complement)

    dvt_direct = dvt[
            ((dvt.allele_1 == dvt.ref) & (dvt.allele_2 == dvt.alt)) |
            ((dvt.allele_1 == dvt.alt) & (dvt.allele_2 == dvt.ref))
            ]

    dvt_complement = dvt[
            ((dvt.allele_1_complement == dvt.ref) & (dvt.allele_2_complement == dvt.alt)) |
            ((dvt.allele_1_complement == dvt.alt) & (dvt.allele_2_complement == dvt.ref))
            ]

    print(f"Strand flipped {len(dvt_complement)} SNPs")

    assert len(
        set(dvt_direct.chrom.astype(str) + "_" + dvt_direct.pos.astype(str)
            ).intersection(
        set(dvt_complement.chrom.astype(str) + "_" + dvt_complement.pos.astype(str)))
        ) == 0

    dvt_complement["allele_1"] = dvt_complement["allele_1_complement"]
    dvt_complement["allele_2"] = dvt_complement["allele_2_complement"]

    dvt = pd.concat([dvt_direct, dvt_complement])
    dvt.drop(columns=['allele_1_complement', 'allele_2_complement'], inplace=True)

    n_snps_missing = len(set(dt.chrom.astype(str) + "_" + dt.pos.astype(str))
                    ) - len(set(dvt.chrom.astype(str) + "_" + dvt_direct.pos.astype(str)))

    if n_snps_missing > 0:
        print(f"Missing {n_snps_missing} SNPs")
    else:
        print("No SNPs missing")
        

    ensure_dir("data/1kg_arrays")
    save_df(dvt, f"1kg_arrays/{array_name}")


@cache
def load_1kg_array(array_name):
    dvt = load_df(f"1kg_arrays/{array_name}")
    assert all(((dvt.allele_1 == dvt.ref) & (dvt.allele_2 == dvt.alt) | 
        (dvt.allele_1 == dvt.alt) & (dvt.allele_2 == dvt.ref)))
    return dvt


def store_ukb_pca_array():
    dt = load_ukb_pca_map()
    store_1kg_array(dt, "ukb_pca")


def get_complement(allele):
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return complement.get(allele, allele)



def map_ancestral_chunk_to_member_1kg(ancestor, ancestor_member):
    name_to_id_1kg = load_name_to_id_1kg()

    suffix = 1 if ancestor_member == "paternal" else 2
    id_1kg = name_to_id_1kg[ancestor]
    member_1kg = f"{id_1kg}_{suffix}"
    return member_1kg


def load_person_dvt(generation, name, array_name):
    dvt = load_1kg_array(array_name)
    
    dfa = load_ancestry(generation=generation)
    dfa = dfa[dfa.name == name]
    
    obligatory_cols = ['chrom', 'pos', 'allele_1', 'allele_2', 'ref', 'alt']
    
    data_paternal = []
    data_maternal = []
    for _, chunk in dfa.iterrows():

        member_1kg = map_ancestral_chunk_to_member_1kg(chunk.ancestor, chunk.ancestor_member)
        
        overlap_mask = (dvt.chrom == chunk.chrom) & (dvt.pos.between(chunk.pos_start, chunk.pos_end))
        
        dvt_chunk = dvt.loc[overlap_mask, obligatory_cols + [member_1kg]]
        dvt_chunk = dvt_chunk.rename(columns={member_1kg: chunk.member})
        
        if chunk.member == "paternal":
            data_paternal.append(dvt_chunk)
        else:
            data_maternal.append(dvt_chunk)

    dvt_paternal = pd.concat(data_paternal, ignore_index=True)
    dvt_maternal = pd.concat(data_maternal, ignore_index=True)

    dvt = dvt_paternal.merge(dvt_maternal, on=obligatory_cols, how='outer')

    return dvt


def get_person_frac_hom(generation, name):
    dvt = load_person_dvt(generation, name, "ukb_pca")
    frac_hom = sum((dvt.paternal == dvt.maternal)) / len(dvt)
    return frac_hom


def get_generation_frac_hom(generation):
    names = load_generation_names(generation)
    data = []
    for name in names:
        frac_hom = get_person_frac_hom(generation, name)
        entry = {
            "generation": generation,
            "name": name,
            "frac_hom": frac_hom
        }
        data.append(entry)
    return pd.DataFrame(data)


def store_homozygosity(ncores=9):
    tasks = list(range(N_GENERATIONS + 1))

    with mp.Pool(processes=ncores) as pool:
        data_all = pool.map(get_generation_frac_hom, tasks)

    dh = pd.concat(data_all, ignore_index=True)
    save_df(dh, "frac_hom")














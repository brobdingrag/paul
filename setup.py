from paul import *





def store_initial_generation():
    generation = 0
    
    dfs = load_df("european_samples")
    dfs = dfs[~dfs.related]
    rename_map = {name: i+1 for i, name in enumerate(dfs.name)}

    df = load_df("height_states")

    ensure_dir(f"data/score_genomes/{generation}")
    for name_orig in dfs.name:
        name = rename_map[name_orig]
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







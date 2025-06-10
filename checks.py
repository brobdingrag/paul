from paul import *



def run_checks():
    check_inter_gen_same()
    check_1kg_scores_match()





def check_1kg_scores_match():
    dfs = load_df("european_samples")
    dfs = dfs[~dfs.related]
    rename_map = {name: i+1 for i, name in enumerate(dfs.name)}

    dh = load_df("height_1kg_scores")
    for name_1kg, name in rename_map.items():
        assert np.allclose(score_person(0, name), dh.query("name == @name_1kg").score)

        


def check_inter_gen_same():
    """Test that the two datasets have the same individuals for each generation."""
    dg = load_generation_scores()
    di = load_intergenerational_scores()
    for generation in range(1, N_GENERATIONS + 1):
        names_inter = di[(di.generation == generation - 1) & (di.selected)].name
        names_gen = dg[dg.generation == generation].name
        assert set(names_inter) == set(names_gen)







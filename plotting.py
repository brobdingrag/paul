from paul import *



def plot_parameters():
    plot_male_female_meiosis_maps()
    plot_male_female_n_crossovers()
    plot_pmf_n_embryos()
    

def plot_simulation_results():
    plot_pgs_dist_generations()
    plot_population_size_generations()
    plot_embryo_dists()
    plot_embryo_ranks()


def plot_male_female_meiosis_maps():

    dm = load_df("meioses")
    dm = dm[dm.chrom != "X"]
    dm.chrom = dm.chrom.astype(int)

    fig, axs = golden_fig(ncols=2, scale=1.3, orient='v')
    dm['chrom_j'] = dm.chrom + np.random.uniform(low=-0.3, high=0.3, size=len(dm))
    def get_ax(dm, sex_parent, ax):
        dot_size = 0.01
        parent = "Dad" if sex_parent == "male" else "Mom"
        df = dm[dm.parent == parent]
        color = MALE if sex_parent == "male" else FEMALE
        ax.make_golden(orient='v')
        sns.scatterplot(
            data=df[~df.start],
            x="pos_start_hg38",
            y='chrom_j',
            color=color,
            edgecolor='none',
            s=dot_size,
            ax=ax
        )
        sns.scatterplot(
            data=df[~df.end],
            x="pos_end_hg38",
            y='chrom_j',
            color=color,
            edgecolor='none',
            s=dot_size,
            ax=ax
        )
        ax.set_xticks([0, 100_000_000, 200_000_000], 
                    ['0', '100M', '200M'])
        
        ax.set_title(f"{sex_parent.title()} crossover map", pad=10)
        ax.set_xlabel("Position (bases)")
        ax.set_ylabel("Chromosome")
        ax.set_yticks(list(range(1, 23)))
        ax.set_ylim(0.5, 22.5)
        ax.reverse_ylim()
        ax.remove_top_right_spines()
        # ax.set_xlim(-5_000_000, 261315433)
        ax.set_xlim(0, 261315433)

    axs[0] = get_ax(dm, "male", axs[0])
    axs[1] = get_ax(dm, "female", axs[1])
    save_fig(fig, "male_female_meioses", pad_inches=0.5, tight_layout=True, w_pad=4)



def plot_male_female_n_crossovers():
    dm = load_df("meioses")
    n_male = ((dm[(dm.parent == "Dad") & (dm.chrom != "X")].groupby("individual").size()).mean() - 22) * 23 / 22
    n_female = ((dm[(dm.parent == "Mom") & (dm.chrom != "X")].groupby("individual").size()).mean() - 22) * 23 / 22
    n_male = np.round(n_male, 2)
    n_female = np.round(n_female, 2)
    n_male_standard = 27
    n_female_standard = 43
    print(f"Male meioses {n_male} vs {n_male_standard}")
    print(f"Female meioses {n_female} vs {n_female_standard}")

    n_male_events = dm[(dm.parent == "Dad") & (dm.chrom != "X")].groupby("individual").size() - 22
    n_female_events = dm[(dm.parent == "Mom") & (dm.chrom != "X")].groupby("individual").size() - 22
    fig, ax = square_fig()
    ax.grid(axis='y')
    sns.histplot(n_male_events, ax=ax, color=MALE, binwidth=1, alpha=0.6, zorder=10)
    sns.histplot(n_female_events, ax=ax, color=FEMALE, binwidth=1, alpha=0.6, zorder=10)
    ax.set_ylabel("Count")
    ax.set_xlim(15, 65)
    ax.remove_top_right_spines()
    ax.comma_yscale()
    ax.set_xlabel("Number of crossovers per meiosis")
    save_fig(fig, "meioses_events", pad_inches=0.25, tight_layout=True)


def load_palette(n_generations=N_GENERATIONS):
    palette = sns.color_palette("viridis", n_colors=n_generations+3)[:n_generations+1] 
    return palette


def plot_population_size_generations():
    sizes = load_population_sizes()
    palette = load_palette()
    
    fig, ax = square_fig()

    ax.plot(range(N_GENERATIONS + 1), sizes, color="black", alpha=0.2, zorder=0)

    def _get_ha(i):
        if i == 0:
            return "center"
        elif i == 1:
            return "center"
        else:
            return "right"

    for i in range(N_GENERATIONS + 1):
        ax.scatter(x=i, y=sizes[i], color=palette[i], zorder=10)
        ax.text(x=i, y=sizes[i] + 130, s=comma(sizes[i]), 
                ha=_get_ha(i), va="bottom", zorder=10,
                fontsize=10, color=palette[i])

    # Design
    ax.remove_top_right_spines()
    ax.grid(axis='y')

    # Y axis
    ax.set_ylabel("Population size")
    ax.comma_yscale()
    ax.set_ylim(0.1, max(sizes) + 500)

    # X axis
    ax.set_xlabel("Generations of selection")
    edge_x = 0.4
    ax.set_xlim(-edge_x, N_GENERATIONS + edge_x)

    save_fig(fig, "population_size_generations", pad_inches=0.25)


def plot_pgs_dist_generations(n_generations=N_GENERATIONS):
    dg = load_generation_scores()
    loc, scale = get_mad_normal(dg[dg.generation == 0].score)
    dg.score = (dg.score - loc) / scale

    fig, ax = golden_fig()
    palette = load_palette(n_generations)
    lw=1.2
    sns.kdeplot(
        data=dg[dg.generation <= n_generations],
        x="score",
        hue="generation",
        common_grid=True,
        common_norm=False,
        bw_adjust=2,
        palette=palette,
        linewidth=lw,
        alpha=0.75,
        fill=False,
    )

    # Manually create the legend for the generations

    # Create custom legend handles
    def get_label(generation):
        if generation == 0:
            return f"0 (1kg EUR)"
        return generation

    handles = [
        Line2D([0], [0], color=palette[i], lw=lw, label=get_label(i))
        for i in range(n_generations+1)
    ]

    ax.legend(handles=handles, title="Generations\nof selection", frameon=False, loc="upper left") 

    ax.naked_top()
    ax.set_xticks(range(-6, 10, 2))
    ax.set_xlim(-7, 9)
    ax.set_ylim(0.001, None)
    ax.sign_xscale()
    ax.set_xlabel("Polygenic score (standard deviations in 1kg EUR)")
    save_fig(fig, f"pgs_dist_generations", pad_inches=0.25)


def plot_pmf_n_embryos():
    dp = load_pmf_n_embryos()

    fig, ax = square_fig()

    sns.barplot(x="n_embryos", y="proba", data=dp, ax=ax, 
                width=1, edgecolor='black', linewidth=0.5,
                color='lightgrey', zorder=10)
    
    # Design
    ax.grid(axis='y')
    ax.remove_top_right_spines()

    # Y axis
    ax.set_ylabel("Probability")
    ax.percent_yscale(change_lim=False)
    ax.set_ylim(0, 0.12)

    # X axis
    ax.set_xlabel("Number of viable embryos per couple")
    xlim = 26
    ax.set_xticks(range(0, 30, 5))
    edge_x = 0.5
    ax.set_xlim(-edge_x, xlim+edge_x)

    save_fig(fig, "pmf_n_embryos", pad_inches=0.25)


def plot_embryo_dists():
    dh = load_df("height_1kg_scores")
    dh = dh[dh.trio.notna()]
    dh.trio = dh.trio.astype(int)
    parental_midpoints = dh[dh.parent].groupby('trio').score.mean()
    child_scores = dh[dh.child].groupby('trio').score.first()
    child_adjusted_scores = child_scores - parental_midpoints 

    di = load_intergenerational_scores()
    dg = load_generation_scores()

    _, scale = get_mad_normal(dg[dg.generation == 0].score)

    fig, axs = square_fig(nrows=3, ncols=3)

    def _get_title(generation):
        if generation == -1:
            return "1kg EUR actual trios (106)"
        else:
            return f"Generation {generation} trios"

    for generation, ax in zip(range(-1, N_GENERATIONS), axs):

        if generation == -1:
            scores_embryos = child_adjusted_scores / scale

        else:
            scores_embryos = di[di.generation == generation].embryo_score_normalized / scale

        binrange = (-2.5, 2.5)
        binwidth = 0.1
        color_main = 'lightgrey' if generation >= 0 else 'purple'
        sns.histplot(scores_embryos, ax=ax, color=color_main, 
                     stat='count',
                     binrange=binrange, binwidth=binwidth)

        legend_options = {"borderpad": 0.6, "bbox_to_anchor": (0, 0.95), "loc": "upper left"}
        if generation == -1:
            ax.legend(["Actual\nchildren"], **legend_options)
        else:
            scores_babies = di[(di.generation == generation) & (di.selected)].embryo_score_normalized / scale
            sns.histplot(scores_babies, ax=ax, color="purple", 
                         stat='count',
                         binrange=binrange, binwidth=binwidth)
            ax.legend(["Embryos", "Babies"], **legend_options)

        # X axis
        ax.set_xlim(binrange)
        ax.set_xlabel("PGS (Cent.)")
        ax.sign_xscale()

        # Y axis
        ax.comma_yscale()
        ax.set_nyticks(5)

        # Design
        # ax.naked_top()
        ax.remove_top_right_spines()
        ax.set_title(_get_title(generation))

        ax.text(0.04, 0.72, s=f"Std. dev. {scores_embryos.std():.2f}",
                transform=ax.transAxes,
                fontsize=10, color="black", alpha=0.5)

    save_fig(fig, "embryo_dists", pad_inches=0.4, h_pad=4)



def plot_embryo_ranks():
    di = load_intergenerational_scores()
    counts = np.zeros(30)
    for generation in range(N_GENERATIONS):
        dt = di[di.generation == generation]
        for couple_num in range(1, dt.couple_num.max() + 1):
            dtc = dt[dt.couple_num == couple_num]
            dtc = dtc.sort_values("embryo_score", ascending=False).reset_index()
            ranks = dtc[dtc.selected].index.tolist()
            for rank in ranks:
                counts[rank] += 1
    fracs = counts / sum(counts)
    dc = pd.DataFrame({"embryo_rank": list(range(1, 31)), "births_frac": fracs})

    fig, ax = square_fig()
    sns.barplot(x="embryo_rank", y="births_frac", data=dc, ax=ax, 
                width=1, edgecolor='black', linewidth=0.5,
                color='lightgrey', zorder=10)
        
    ax.set_xlim(-0.5, 9.5)
    ax.set_xticks(range(0, 10))
    # ax.set_xlim(-0.5, 10.5)

    ax.percent_yscale(change_lim=False)
    ax.set_xlabel("Embryo rank")
    ax.set_ylabel("Fraction of births")

    ax.grid(axis='y')
    ax.set_ylim(0, 0.24)
    ax.remove_top_right_spines()
    save_fig(fig, "embryo_ranks", pad_inches=0.25)




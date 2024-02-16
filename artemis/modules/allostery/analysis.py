import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import itertools as itt
from scipy.stats import zscore
from scipy.special import comb
from tqdm import tqdm

width = 2
top = None
width2 = .6
max_iters = 10000

print("\nSCRIPT FOR ALLOSTERY ANALYSIS IS LAUNCHED\n")

def analyze_allostery(Allostery, out_path, top, noseq, Zscore):

    if Allostery.allosteric_site is None:
        print("It is not possible to perform the analysis without an allosteric site. Interrupted.\n")

    intensity = []
    for i in range(Allostery.NResidues):
        if i + 1 in Allostery.active_site:
            intensity.append(-1)
        else:
            inten = 0
            for resid in Allostery.active_site:
                if np.abs(resid - 1 - i) >= noseq:
                    inten += Allostery.map_[resid - 1][i]
            intensity.append(inten)

    NAct = len(Allostery.active_site)
    NAll = len(Allostery.allosteric_site)

    act_intensity = np.zeros(NAct)
    all_intensity = np.zeros(NAll)
    act_entropy = Allostery.map_[Allostery.active_site-1, Allostery.active_site-1]
    all_entropy = Allostery.map_[Allostery.allosteric_site-1, Allostery.allosteric_site-1]
    act_Entropy = 0
    all_Entropy = 0
    connectivity = 0
    act_names = []
    all_names = []
    for i, res in enumerate(Allostery.active_site):
        act_names.append(Allostery.names[res-1] + "\n(" + str(Allostery.real_numbers[res-1]) +")")
        for j in range(Allostery.NResidues):
            if not (j+1 in Allostery.active_site):
                act_intensity[i] += Allostery.map_[res-1, j]
            elif (j+1 == res):
                act_Entropy += Allostery.map_[res-1, j]
            else:
                act_Entropy += 1/2*Allostery.map_[res-1, j]

            if (j+1 in Allostery.allosteric_site):
                connectivity += Allostery.map_[res-1, j]

    for i, res in enumerate(Allostery.allosteric_site):
        all_names.append(Allostery.names[res-1] + "\n(" + str(Allostery.real_numbers[res-1]) +")")
        for j in range(Allostery.NResidues):
            if not (j+1 in Allostery.allosteric_site):
                all_intensity[i] += Allostery.map_[res-1, j]
            elif (j+1 == res):
                all_Entropy += Allostery.map_[res-1, j]
            else:
                all_Entropy += 1/2*Allostery.map_[res-1, j]

    act_mi = np.sum(act_intensity)
    all_mi = np.sum(all_intensity)

    size = np.max([NAct, NAll])

    print('\n')

    print("The mutual information of the active site with the remaining protein is", round(act_mi, 2))
    print("The mutual information of the allosteric site with the remaining protein is", round(all_mi, 2))
    print("The internal two-dimensional entropy of the active site is", round(act_Entropy, 2))
    print("The internal two-dimensional entropy of the allosteric site is", round(all_Entropy, 2))
    print("Allosteric connectivity between active and allosteric sites is", round(connectivity, 2), '\n\n')

    fig = plt.figure(figsize=(size*width + 5, 24))
    axs = [None for _ in range(4)]

    axs[0] = plt.subplot2grid((2,2), (0,0),colspan=1, rowspan=1)
    INTENSITY = {"Intensity": np.array(act_intensity), "Residue": act_names}
    colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
    colors = (sns.color_palette(colors, as_cmap=True)).tolist()
    axs[0] = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, ax=axs[0], hue="Residue", legend=False)
    axs[0].set_title('Intensity of connectivity of residues\nfrom the active site', fontsize=30)
    axs[0].tick_params(axis='both', which='major', labelsize=16)
    axs[0].set_ylabel("Intensity", fontsize=20)

    axs[1] = plt.subplot2grid((2,2), (0,1),colspan=1, rowspan=1)
    INTENSITY = {"Intensity": np.array(all_intensity), "Residue": all_names}
    colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
    colors = (sns.color_palette(colors, as_cmap=True)).tolist()
    axs[1] = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, ax=axs[1], hue="Residue", legend=False)
    axs[1].set_title('Intensity of connectivity of residues\nfrom the allosteric site', fontsize=30)
    axs[1].tick_params(axis='both', which='major', labelsize=16)
    axs[1].set_ylabel("Intensity", fontsize=20)

    axs[2] = plt.subplot2grid((2,2), (1,0),colspan=1, rowspan=1)
    INTENSITY = {"Intensity": np.array(act_entropy), "Residue": act_names}
    colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
    colors = (sns.color_palette(colors, as_cmap=True)).tolist()
    axs[2] = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, ax=axs[2], hue="Residue", legend=False)
    axs[2].set_title('Two-dimensional entropy of residues\nfrom the active site', fontsize=30)
    axs[2].tick_params(axis='both', which='major', labelsize=16)
    axs[2].set_ylabel("Intensity", fontsize=20)

    axs[3] = plt.subplot2grid((2,2), (1,1),colspan=1, rowspan=1)
    INTENSITY = {"Intensity": np.array(all_entropy), "Residue": all_names}
    colors = plt.cm.viridis(INTENSITY["Intensity"]/np.max(INTENSITY["Intensity"]))
    colors = (sns.color_palette(colors, as_cmap=True)).tolist()
    axs[3] = sns.barplot(x="Residue", y="Intensity", data=INTENSITY, palette=colors, dodge=False, ax=axs[3], hue="Residue", legend=False)
    axs[3].set_title('Two-dimensional entropy of residues\nfrom the allosteric site', fontsize=30)
    axs[3].tick_params(axis='both', which='major', labelsize=16)
    axs[3].set_ylabel("Intensity", fontsize=20)

    try:
        fig.savefig(out_path)
        print("File",out_path + " created")
    except:
        print("Error writing file",out_path)


    if Zscore:
        intensity = []
        Allosteric = []
        for i in range(Allostery.NResidues):

            if i+1 in Allostery.allosteric_site:
                Allosteric.append(True)
            else:
                Allosteric.append(False)

            if i + 1 in Allostery.active_site:
                intensity.append(0)
            else:
                inten = 0
                for resid in Allostery.active_site:
                    if np.abs(resid - 1 - i) >= noseq:
                        inten += Allostery.map_[resid - 1][i]
                intensity.append(inten)

        sigma = np.std(intensity)
        mean = np.mean(intensity)

        intensity = zscore(intensity)

        new_names = []
        for i in range(Allostery.NResidues):
            new_names.append(Allostery.names[i] + " (" + str(Allostery.real_numbers[i]) +")")

        unique, counts = np.unique(new_names, return_counts=True)

        for i, item in enumerate(unique):
            if counts[i] > 1:
                indeces = np.argwhere(np.array(new_names)==item)
                for j, idx in enumerate(indeces):
                    new_names[int(idx)] += "[" +str(j+1)+"]"

        INTENSITY = {}
        INTENSITY["Intensity"] = np.array(intensity)
        INTENSITY["Residue"] = new_names
        INTENSITY["Allosteric site"] = Allosteric
        INTENSITY = pd.DataFrame(data=INTENSITY)
        INTENSITY = INTENSITY.sort_values(by=['Intensity'], ascending=False, ignore_index=True)
        try:
            fname = ''.join(out_path.split('.')[:-1]) + '_zscore.csv'
            INTENSITY.to_csv(fname, index=False)
            print("\nFile",fname + " created\n")
        except:
            print("\nError writing file",fname +'\n')


        NSites = len(Allostery.active_site)

        cut = NSites
        iters = 0
        for i in range(NSites):
            if iters >= max_iters:
                cut = i
                print("Warning (zscore_top.py): The size of the active site is too big. Combinations of at most " + str(cut) + " remainders are investigated.\n")
                break
            iters += comb(NSites, i+1)

        Combinations = []
        for i in range(NSites):
            combinations = itt.combinations(Allostery.active_site, i+1)
            Combinations.append(combinations)


        table = pd.DataFrame(columns=['Combination', 'Residues with zscore > 2', 'Zscore > 2 Residues names', 'Residues with zscore > 1.5',
                                    'Zscore > 1.5 Residues names', 'Residues with zscore > 1', 'Zscore > 1 Residues names',  'Mean Intensity', 'STD',
                                    'top 1 Residue', 'top 2 Residue', 'top 3 Residue'])

        more2_ = []
        more1_all = []
        columns = ["Residue"]
        leng = 0

        for combinations in tqdm(Combinations, desc ="Progress", total=cut):
            leng += 1
            if leng > cut:
                break

            for active_site in tqdm(combinations, desc =f"ActSite size {leng} progress", total=comb(NSites, leng)):
                intensity = []
                Allosteric = []
                for i in range(Allostery.NResidues):

                    if i+1 in Allostery.allosteric_site:
                        Allosteric.append(True)
                    else:
                        Allosteric.append(False)

                    if i + 1 in Allostery.active_site:
                        intensity.append(0)
                    else:
                        inten = 0
                        for resid in Allostery.active_site:
                            if np.abs(resid - 1 - i) >= noseq:
                                inten += Allostery.map_[resid - 1][i]
                        intensity.append(inten)

                sigma = np.std(intensity)
                mean = np.mean(intensity)

                intensity = zscore(intensity)

                new_names = []
                for i in range(Allostery.NResidues):
                    new_names.append(Allostery.names[i] + " (" + str(Allostery.real_numbers[i]) +")")

                unique, counts = np.unique(new_names, return_counts=True)

                for i, item in enumerate(unique):
                    if counts[i] > 1:
                        indeces = np.argwhere(np.array(new_names)==item)
                        for j, idx in enumerate(indeces):
                            new_names[int(idx)] += "[" +str(j+1)+"]"

                INTENSITY = {}
                INTENSITY["Intensity"] = np.array(intensity)
                INTENSITY["Residue"] = new_names
                INTENSITY["Allosteric site"] = Allosteric
                INTENSITY = pd.DataFrame(data=INTENSITY)
                INTENSITY = INTENSITY.sort_values(by=['Intensity'], ascending=False, ignore_index=True)

                top1 = str(INTENSITY['Residue'][0]) + " [" +  str(np.round(INTENSITY['Intensity'][0], 2)) + "]"
                top2 = str(INTENSITY['Residue'][1]) + " [" +  str(np.round(INTENSITY['Intensity'][1], 2)) + "]"
                top3 = str(INTENSITY['Residue'][2]) + " [" +  str(np.round(INTENSITY['Intensity'][2], 2)) + "]"

                act_s_name = []
                for i in Allostery.active_site:
                    act_s_name.append(Allostery.names[i-1] + " (" + str(Allostery.real_numbers[i-1]) +")")

                unique, counts = np.unique(act_s_name, return_counts=True)

                for i, item in enumerate(unique):
                    if counts[i] > 1:
                        indeces = np.argwhere(np.array(act_s_name)==item)
                        for j, idx in enumerate(indeces):
                            act_s_name[int(idx)] += "[" +str(j+1)+"]"

                m2R = (INTENSITY[INTENSITY["Intensity"] > 2]["Allosteric site"] == True).sum(axis=0)
                m15R = (INTENSITY[INTENSITY["Intensity"] > 1.5]["Allosteric site"] == True).sum(axis=0)
                m1R = (INTENSITY[INTENSITY["Intensity"] > 1]["Allosteric site"] == True).sum(axis=0)
                INTENSITY1 = INTENSITY[INTENSITY["Intensity"] > 1]
                INTENSITY15 = INTENSITY[INTENSITY["Intensity"] > 1.5]
                INTENSITY2 = INTENSITY[INTENSITY["Intensity"] > 2]
                act_s1 = np.array(INTENSITY1[INTENSITY1["Allosteric site"] == True]["Residue"])
                act_s15 = np.array(INTENSITY15[INTENSITY15["Allosteric site"] == True]["Residue"])
                act_s2 = np.array(INTENSITY2[INTENSITY2["Allosteric site"] == True]["Residue"])

                new_row = pd.Series({'Combination': ", ".join(map(str,act_s_name)), 'Residues with zscore > 2': m2R, 'Zscore > 2 Residues names': ", ".join(map(str,act_s2)),
                                    'Residues with zscore > 1.5': m15R, 'Zscore > 1.5 Residues names': ", ".join(map(str,act_s15)),'Residues with zscore > 1': m1R,
                                    'Zscore > 1 Residues names': ", ".join(map(str,act_s1)), 'Mean Intensity': np.round(mean, 2), 'STD': np.round(sigma, 2),
                                    'top 1 Residue': top1, 'top 2 Residue': top2, 'top 3 Residue': top3})

                table = pd.concat([table, new_row.to_frame().T], ignore_index=True)

                if len(active_site) == 1:
                    more2_.append(np.array(INTENSITY2["Residue"]))
                    more1_all.append(np.array(INTENSITY1[INTENSITY1["Allosteric site"] == True]["Residue"]))
                    columns.append(", ".join(map(str,act_s_name)))
            print('\033[F\033[K\033[F')

        try:
            fname = ''.join(out_path.split('.')[:-1]) + '_zscore_top.csv'
            table.to_csv(fname, index=False)
            print("\nFile",fname + " created\n")
        except:
            print("\nError writing file",fname + '\n')

        overlap = pd.DataFrame(columns=columns)
        overlap_all = pd.DataFrame(columns=columns)
        ovAmounts = pd.DataFrame(columns=columns)
        ovAmounts_all = pd.DataFrame(columns=columns)
        columns = columns[1:]
        overlap["Residue"] = columns
        overlap_all["Residue"] = columns
        ovAmounts["Residue"] = columns
        ovAmounts_all["Residue"] = columns

        for resid1 in range(NSites):
            lists = []
            lists_all = []
            amounts = []
            amounts_all = []
            for resid2 in range(NSites):
                list_ = list(set(more2_[resid1]) & set(more2_[resid2]))
                list_all = list(set(more1_all[resid1]) & set(more1_all[resid2]))
                lists.append(", ".join(map(str,list_)))
                lists_all.append(", ".join(map(str,list_all)))
                amounts.append(len(np.array(list_)))
                amounts_all.append(len(np.array(list_all)))

            overlap[columns[resid1]] = lists
            overlap_all[columns[resid1]] = lists_all
            ovAmounts[columns[resid1]] = amounts
            ovAmounts_all[columns[resid1]] = amounts_all

        try:
            fname = ''.join(out_path.split('.')[:-1]) + '_overlap.xlsx'
            with pd.ExcelWriter(fname) as writer:
                overlap.to_excel(writer, index=False, sheet_name='overlap (lists)')
                overlap_all.to_excel(writer, index=False, sheet_name='Allosteric sites overlap (lists)')
                ovAmounts.to_excel(writer, index=False, sheet_name='overlap (amount)')
                ovAmounts_all.to_excel(writer, index=False, sheet_name='Allosteric sites overlap (amount)')
            print("\nFile",fname + " created\n")
        except:
            print("\nError writing file",fname + '\n')


    if top is not None:
        NSites = len(Allostery.active_site)

        cut = NSites
        iters = 0
        for i in range(NSites):
            if iters >= max_iters:
                cut = i
                print("Warning (intensity_top.py): The size of the active site is too big. Combinations of at most " + str(cut) + " remainders are investigated.\n")
                break
            iters += comb(NSites, i+1)

        Combinations = []
        for i in range(NSites):
            combinations = itt.combinations(Allostery.active_site, i+1)
            Combinations.append(combinations)


        table = pd.DataFrame(columns=['Combination', 'Top Residues', 'Top Residues names', 'Mean Intensity', 'STD',
                                    'top 1 Residue', 'top 2 Residue', 'top 3 Residue'])

        leng = 0
        for combinations in tqdm(Combinations, desc ="Progress", total=cut):
            leng += 1
            if leng > cut:
                break

            for active_site in tqdm(combinations, desc =f"ActSite size {leng} progress", total=comb(NSites, leng)):
                intensity = []
                Allosteric = []
                for i in range(NResidues):

                    if i+1 in Allostery.allosteric_site:
                        Allosteric.append(True)
                    else:
                        Allosteric.append(False)

                    if i + 1 in Allostery.active_site:
                        intensity.append(0)
                    else:
                        inten = 0
                        for resid in Allostery.active_site:
                            if np.abs(resid - 1 - i) >= noseq:
                                inten += Allostery.map_[resid - 1][i]
                        intensity.append(inten)

                new_names = []
                for i in range(Allostery.NResidues):
                    new_names.append(Allostery.names[i] + " (" + str(Allostery.real_numbers[i]) +")")

                unique, counts = np.unique(new_names, return_counts=True)

                for i, item in enumerate(unique):
                    if counts[i] > 1:
                        indeces = np.argwhere(np.array(new_names)==item)
                        for j, idx in enumerate(indeces):
                            new_names[int(idx)] += "[" +str(j+1)+"]"

                INTENSITY = {}
                INTENSITY["Intensity"] = np.array(intensity)
                INTENSITY["Residue"] = new_names
                INTENSITY["Allosteric site"] = Allosteric
                INTENSITY = pd.DataFrame(data=INTENSITY)
                INTENSITY = INTENSITY.sort_values(by=['Intensity'], ascending=False, ignore_index=True)

                top10 = int(top/100*NResidues)

                top1 = str(INTENSITY['Residue'][0]) + " [" +  str(np.round(INTENSITY['Intensity'][0], 2)) + "]"
                top2 = str(INTENSITY['Residue'][1]) + " [" +  str(np.round(INTENSITY['Intensity'][1], 2)) + "]"
                top3 = str(INTENSITY['Residue'][2]) + " [" +  str(np.round(INTENSITY['Intensity'][2], 2)) + "]"

                act_s_name = []
                for i in active_site:
                    act_s_name.append(Allostery.names[i-1] + " (" + str(Allostery.real_numbers[i-1]) +")")

                unique, counts = np.unique(act_s_name, return_counts=True)

                for i, item in enumerate(unique):
                    if counts[i] > 1:
                        indeces = np.argwhere(np.array(act_s_name)==item)
                        for j, idx in enumerate(indeces):
                            act_s_name[int(idx)] += "[" +str(j+1)+"]"

                t10R = (INTENSITY.iloc[:top10]["Allosteric site"] == True).sum(axis=0)
                INTENSITY = INTENSITY.iloc[:top10]
                act_s = np.array(INTENSITY[INTENSITY["Allosteric site"] == True]["Residue"])

                new_row = pd.Series({'Combination': ", ".join(map(str,act_s_name)), 'Top Residues': t10R, 'Top Residues names': ", ".join(map(str,act_s)),
                                    'Mean Intensity': np.round(np.mean(INTENSITY["Intensity"]),2), 'STD': np.round(np.std(INTENSITY["Intensity"]), 2),
                                    'top 1 Residue': top1, 'top 2 Residue': top2, 'top 3 Residue': top3})

                table = pd.concat([table, new_row.to_frame().T], ignore_index=True)
            print('\033[F\033[K\033[F')

        try:
            fname = ''.join(out_path.split('.')[:-1]) + "_top"+str(top)+".csv"
            table.to_csv(fname, index=False)
            print("\nFile",fname+" created\n")
        except:
            print("\nError writing file",fname+'\n')

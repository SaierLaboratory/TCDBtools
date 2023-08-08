#!/usr/bin/env python

def figHydroAlignedRegions(seq,seqfrag,sbj,sbjfrag,seqname,sbjname,savefig):
    from Bio import Seq, SeqIO
    import matplotlib.pyplot as plt
    import libquod
    
    fig, ax = plt.subplots()
    
    #seq1 = SeqIO.SeqRecord(seq=Seq.Seq("MNKQKNVSKHPPGLYLLFFTEMWERFSYYGMRGILVLYLTKELIEGGLGMSSTVALSIYGFYTGACYLAPLVGGWLSDMYLGKRLAITLGGISMALGNFIMFLGDAQWSMYTGLALLIIGNGFFKPNISTLVGELYEENDPKRDAAFTIFYMGINLGAFIAPLIIALVTDKIMGGMQGDVFVYGYKYGFLISSIGMVIGQVLFNTLGNKYLGDLGKKPSRNEVVGDVEAKNRPLTSIEKQRTVAILILASIVIFFWAGFEQAGSSLTLYTDKFVDRTIFGFEIPTPWFQSVNPLFIILLAPLVSMLWVKLARSKRGDLPSPVKMAIGMILLGVGYMVLLVAVLGTGSNEAEVVNKAHIIFIVLTYMFHTIGELFLSPVGLSTVSRIAPVKLASLLMGVWLASSGIANILAGQLASVTQKLGYFEVFAAIGGMAILLGLILLAFSKKITKMMHID"), id="WP_019240224", name="", description="")
    #fragseq1 = SeqIO.SeqRecord(seq=Seq.Seq("IEGGLGMSSTVALSIYGFYTGACYLAPLVGGWLSDMYLGKRLAITLG-GISMALGNFIMFLGDAQWSMYTGLALLI-IGNGFFKPNISTLVGELYEENDPKRDAAFTIFYMGINLGAFIAPLIIALVTDKIMGGMQGDVFVYGYKYGFLISSI-GMVIGQVLFNTLGNKYLGDLGKKPSRNEVVGDVEAKNRPLTSIEKQRTVAILILASIVIFFWAGFEQAGSSLTLYTDKF-VDRTIFGFEIPTPWFQSVNPLFIILLAPLVSMLWVKLARSKRGDLPSPVKMAIGMILLGVGYMVLLVAVLGTGSNEAEVVNKAHIIFIVLTYMFHTIGELFLSPVGLSTVSRIAPVKLASLLMGVWLASSGIANILAGQLASVTQKLGYFEVFAAIGGMAILLGLILLA"), id="WP_019240224_frag", name="", description="")
    
    seq1 = SeqIO.SeqRecord(seq=Seq.Seq(seq), id=f"{seqname}", name=f"{seqname}", description="")
    fragseq1 = SeqIO.SeqRecord(seq=Seq.Seq(seqfrag), id=f"{seqname}_frag", name="", description="")
    
    seq2 = SeqIO.SeqRecord(seq=Seq.Seq(sbj), id=f"{sbjname}", name="", description="")
    fragseq2 = SeqIO.SeqRecord(seq=Seq.Seq(sbjfrag), id=f"{sbjname}_frag", name="", description="")
    #Set the figure width and heigth
    fig.set_figwidth(10)
    fig.set_figheight(4)
    
    #to get rid of spaces to the left, right and top of the plots
    fig.set_tight_layout(True)
    
    def plot_what(seq, fragseq=None, ax=None, fc=None, ec=None):
        out = {
                "hydro":libquod.entities.Hydropathy.compute(seq, fragment=fragseq, ec=ec),
                "hmmtop":libquod.entities.HMMTOP.compute(seq, fragment=fragseq, fc=fc),
                }
        if ax is None:
            fig, ax = plt.subplots()
    
        out["hydro"].plot(ax=ax)
        out["hmmtop"].plot(ax=ax)
    
        return out
    
    #plot_what(seq, ec="red", fc="orange", ax=ax)
    plot_what(seq1, fragseq1, ec="red", fc="orange", ax=ax)
    plot_what(seq2, fragseq2, ec="blue", fc="cyan", ax=ax)
    
    ax.set_xlim([0, len(fragseq1)])
    ax.set_ylim([-3, 3])
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    #The xticks locations
    xticks = list(range(0, len(fragseq1), 50))
    ax.set_xticks(xticks)
    
    plt.grid()
    
    #The axes and plot labels.
    plt.title(f'{seqname} (red) vs {sbjname} (blue)')
    plt.xlabel('Residues')
    plt.ylabel('Hydropathy (kcal/mol)')
    plt.savefig(savefig)
    
    #plt.show()

def figFullProtein(seq,seqname,savefig,chart3,startend):
    from Bio import Seq, SeqIO
    
    import matplotlib.pyplot as plt
    import libquod
    
    fig, ax = plt.subplots()
    
    seq = SeqIO.SeqRecord(seq=Seq.Seq(seq), id="test",  name=f"{seqname}", description="")
    
    #Set the figure width and heigth
    fig.set_figwidth(10)
    fig.set_figheight(4)
    
    #to get rid of spaces to the left, right and top of the plots
    fig.set_tight_layout(True)
    
    #Calculate ticks in the x-axis
    xticks = list(range(0, len(seq), 50))
    
    
    #single sequence
    libquod.entities.Hydropathy.compute(seq, edgecolor='red').plot(ax=ax)
    libquod.entities.HMMTOP.compute(seq, fc="orange").plot(ax=ax)
    
    #Pfam domain: PF12730, which spans 203-377 and is drawn from
    #y=-2.9 to y=-2.75 in red
    if chart3:
        for ind,acc in enumerate(chart3):
            libquod.entities.Region(spans=[[float(acc["from_1"]),float(acc["to_1"])]], yspan=[-2.9+0.5*ind,-2.75+0.5*ind], text=acc["accession"].split(".")[0], valign="t", halign="c", fc="red", alpha=1.0).plot(ax=ax)
    else:
            libquod.entities.Region(spans=[[float(0),float(0)]], yspan=[-2.9,-2.75], text="", valign="t", halign="c", fc="red", alpha=1.0).plot(ax=ax)

        
    #Vertical walls delimiting the region of a protein involved in an alignment.
    #Draw wall enclosing the x-axis coordinates [195, 427]
    libquod.entities.Wall(spans=[startend], y=3.0, ypos="+-", text="Alignment region").plot(ax=ax)
    #libquod.entities.Wall(spans=[hsp.query_start,hsp.query_end], y=3.0, ypos="+-", text="Alignment region").plot(ax=ax)
    
    #Draw horizontal line for hydropathy = 0.0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    #The limits of the y-axis
    ax.set_ylim([-3,3.5])
    
    #The xticks locations
    ax.set_xticks(xticks)
    
    #Draw the grid for the plot
    plt.grid()
    
    
    #The axes and plot labels.
    plt.title(f'{seqname}')
    plt.xlabel('Residues')
    plt.ylabel('Hydropathy (kcal/mol)')
    plt.savefig(savefig)
    
    #plt.show()

def figFullProteinSbj(seq,seqname,savefig,chart2,startendsbj):
    from Bio import Seq, SeqIO
    
    import matplotlib.pyplot as plt
    import libquod
    
    fig, ax = plt.subplots()
    
    seq = SeqIO.SeqRecord(seq=Seq.Seq(seq), id="test",  name=f"{seqname}", description="")
    
    #Set the figure width and heigth
    fig.set_figwidth(10)
    fig.set_figheight(4)
    
    #to get rid of spaces to the left, right and top of the plots
    fig.set_tight_layout(True)
    
    #Calculate ticks in the x-axis
    xticks = list(range(0, len(seq), 50))
    
    
    #single sequence
    libquod.entities.Hydropathy.compute(seq, edgecolor='blue').plot(ax=ax)
    libquod.entities.HMMTOP.compute(seq, fc="cyan").plot(ax=ax)
    
    #Pfam domain: PF12730, which spans 203-377 and is drawn from
    #y=-2.9 to y=-2.75 in blue
    if chart2:
        for ind,acc in enumerate(chart2):
            libquod.entities.Region(spans=[[float(acc["from_1"]),float(acc["to_1"])]], yspan=[-2.9+0.5*ind,-2.75+0.5*ind], text=acc["accession"].split(".")[0], valign="t", halign="c", fc="blue", alpha=1.0).plot(ax=ax)
    else:
            libquod.entities.Region(spans=[[float(0),float(0)]], yspan=[-2.9,-2.75], text="", valign="t", halign="c", fc="blue", alpha=1.0).plot(ax=ax)

    #Vertical walls delimiting the region of a protein involved in an alignment.
    #Draw wall enclosing the x-axis coordinates [195, 427]
    libquod.entities.Wall(spans=[startendsbj], y=3.0, ypos="+-", text="Alignment region").plot(ax=ax)
    
    #Draw horizontal line for hydropathy = 0.0
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    #The limits of the y-axis
    ax.set_ylim([-3,3.5])
    
    #The xticks locations
    ax.set_xticks(xticks)
    
    #Draw the grid for the plot
    plt.grid()
    
    
    #The axes and plot labels.
    plt.title(f'{seqname}')
    plt.xlabel('Residues')
    plt.ylabel('Hydropathy (kcal/mol)')
    plt.savefig(savefig)

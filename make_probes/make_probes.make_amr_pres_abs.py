from Bio import SeqIO

amr_genes = list(SeqIO.parse("CARD_AMR_v3.0.8_curated_Aug2020.fasta", "fasta"))
plas_reps = list(SeqIO.parse("PlasmidFinder_inclZ66.fasta", "fasta"))

#wanted_genes = ["326__Cat_Phe__catA1__2229", "390__Dfr_Tmt__dfrA7__2390", "338__Sul_Sul__sul1__2283", "391__StrA_AGly__strA.v1__2409", "392__StrB_AGly__strB.v1__2410", "338__Sul_Sul__sul2__2284", "282__MphA_MLS__mphA__2174", "143__TEM_Bla__TEM-1__1667", "79__CTX-M_Bla__CTX-M-15__530", "178__QnrS_Flq__qnrS1__2006", "267__Erm_MLS__ermB.v2__2137", "390__Dfr_Tmt__dfrA5__2387", "390__Dfr_Tmt__dfrA15.v2__2365", "368__TetA_Tet__tetA.v2__2318", "369__TetB_Tet__tetB.v2__2325"]

# KH updated 210222
#JH update 07042024 - added NDM-5, KPC-2, OXA-48, VIM-1, IMP-27, removed ampC1 and OXA134
wanted_genes = ["79__CTX-M_Bla__CTX-M-15__530","117__OXA_Bla__OXA-7__927","135__SHV-OKP-LEN_Bla__SHV-12__1581","143__TEM_Bla__TEM-1D.v1__1799","326__Cat_Phe__catA1__2229","327__CmlA_Phe__cmlA1__2259","390__Dfr_Tmt__dfrA14.v1__2362","390__Dfr_Tmt__dfrA15.v2__2365","390__Dfr_Tmt__dfrA17__2369","390__Dfr_Tmt__dfrA5__2387","390__Dfr_Tmt__dfrA7__2390","264__EreA_MLS__ereA__2108","174__QnrB_Flq__qnrB1.v1__1930","176__QnrD_Flq__qnrD1__2002","178__QnrS_Flq__qnrS1__2006","338__Sul_Sul__sul1__2283","338__Sul_Sul__sul2__2284","368__TetA_Tet__tetA.v2__2318","369__TetB_Tet__tetB.v2__2325","370__TetC_Tet__tetC__2330","371__TetD_Tet__tetD__2331","282__MphA_MLS__mphA__2174","267__Erm_MLS__ermB.v2__2137","390__Dfr_Tmt__dfrA19__2371", "390__Dfr_Tmt__dfrA1.v2__2358", "114__NDM_Bla__NDM-5__894", "104__KPC_Bla__KPC-2__815", "117__OXA_Bla__OXA-48__1031", "151__VIM_Bla__VIM-1__1820", "100__IMP_Bla__IMP-27__757"]

#wanted_genes = ["38__AmpC1_Bla__AmpC1__286","79__CTX-M_Bla__CTX-M-15__530","117__OXA_Bla__OXA-7__927","117__OXA_Bla__OXA-134__1022","135__SHV-OKP-LEN_Bla__SHV-12__1581","143__TEM_Bla__TEM-1D.v1__1799","326__Cat_Phe__catA1__2229","327__CmlA_Phe__cmlA1__2259","264__EreA_MLS__ereA__2108","174__QnrB_Flq__qnrB1.v1__1930","176__QnrD_Flq__qnrD1__2002","178__QnrS_Flq__qnrS1__2006","338__Sul_Sul__sul1__2283","338__Sul_Sul__sul2__2284","368__TetA_Tet__tetA.v2__2318","369__TetB_Tet__tetB.v2__2325","370__TetC_Tet__tetC__2330","371__TetD_Tet__tetD__2331","282__MphA_MLS__mphA__2174","267__Erm_MLS__ermB.v2__2137", "390__Dfr_Tmt__dfrA1.v2__2358"]

#JH updated 07042024 - added IncX1 and p0111
wanted_plas_reps = {"125__FIAHI1_HI1__FIAHI1_1_HI1_AF250878__80": "IncFIAHI1", "122__HI1A___HI1A_1__AF250878__74": "IncHI1A", "101__HI1BR27_R27__HI1BR27_1_R27_AF250878__142": "IncHI1BR27", "72__Y_1__Y_1__K02380__313": "IncY", "fljB_z66":"z66", "129__X3___X3_1__JN247852__89": "IncX3", "93__HI2A_1__HI2A_1__BX664015__186": "IncHI2A", "159__I1_Alpha__I1_1_Alpha_AP005147__238": "IncI1", "78__LM___L___M_1__AF550415__320": "IncL_M", "55__FIBpHCM2_pHCM2__FIBpHCM2_1_pHCM2_AL513384__97": "IncFIB_pHCM2", "59__FIBK_Kpn3__FIBK_1_Kpn3_JN233704__92": "IncFIB_K", "108__N___N_1__AY046276__137": "IncN", "126__X1_1__X1_1_EU370913__85": "IncX1", "53__p0111_1__p0111_1__AP010962__107": "p0111"}

probe_seqs = []

for amr_gene in amr_genes:
    if amr_gene.id in wanted_genes:
        # then select this gene to make a probe from
        gene_name = amr_gene.id.split('__')[2]
        # remove the .vX from the gene name
        if '.v' in gene_name:
            gene_name = gene_name.split('.v')[0]
        if '-' in gene_name:
            gene_name = gene_name.replace('-', '')
        if gene_name == 'dfrA19': #CARD lists the pathogenWatch dfrA18 allele as dfrA19, so changing the name here
            gene_name = 'dfrA18'
        # create the mykrobe probe header
        gene_header = gene_name + '?name=' + gene_name + '&version=1'
        amr_probe_record = amr_gene
        amr_probe_record.id = gene_header
        amr_probe_record.description = ''
        probe_seqs.append(amr_probe_record)
for plas_rep in plas_reps:
    if plas_rep.id in wanted_plas_reps.keys():
        # select this rep and make a probe
        probe_name = wanted_plas_reps[plas_rep.id]
        probe_header = probe_name + '?name=' + probe_name + '&version=1'
        plas_probe_record = plas_rep
        plas_probe_record.id = probe_header
        plas_probe_record.description = ''
        probe_seqs.append(plas_probe_record)

#write out
SeqIO.write(probe_seqs, 'typhi.amr.plas.presAbs.20240407.fa', 'fasta')

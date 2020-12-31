# December 2020 
# Seattle, WA
import pandas as pd 
from tcrdist.repertoire import TCRrep
from tcrdist.tabulate import tabulate
import os
import numpy as np
from tcrdist.paths import path_to_base
import time
import math
import parmap

def get_safe_chunk(search_clones, bulk_clones,target = 10**7):
	"""
	This function help pick a chunk size that prevents excessive memory use,
	With two CPU, 10*7 should keep total overall memory demand below 1GB
	"""
	ideal_divisor = (search_clones * bulk_clones) / target
	if ideal_divisor < 1:
		ideal_chunk_size = search_clones
		print(ideal_chunk_size)
	else:
		ideal_chunk_size = math.ceil((search_clones)/ ideal_divisor)
		print(ideal_chunk_size)
	return ideal_chunk_size

def do_search2(file, df_search, dest, tag):
	sample_name = file.replace('.tsv.tcrdist3.v_max.tsv','')
	tic = time.perf_counter()
	
	# <tr_search> tcrdist.repertoire.TCRrep object for computing distances
	tr_search = TCRrep(cell_df = df_search,
					organism = 'human',
					chains = ['beta'],
					db_file = 'alphabeta_gammadelta_db.tsv',
					compute_distances = False)
	# set cpus according to parameter above
	tr_search.cpus = 1
	df_bulk   = pd.read_csv(os.path.join(path, file), sep = '\t')
	df_bulk = df_bulk[['cdr3_b_aa','v_b_gene','j_b_gene','templates','productive_frequency']].rename(columns = {'templates':'count'})

	tr_bulk = TCRrep(	cell_df = df_bulk,                 
						organism = 'human', 
						chains = ['beta'], 
						db_file = 'alphabeta_gammadelta_db.tsv',
						compute_distances = False)
	
	#lines_per_file.append(tr_bulk.clone_df.shape[0]) 
	
	search_clones = tr_search.clone_df.shape[0]
	bulk_clones   = tr_bulk.clone_df.shape[0]
	# To avoid memory pressure on the system we set a target that tcrdist doesn't do more than 10M comparisons per process
	ideal_chunk_size = get_safe_chunk(tr_search.clone_df.shape[0], tr_bulk.clone_df.shape[0],target = 10**8)
	tr_search.compute_sparse_rect_distances(df = tr_search.clone_df, df2 = tr_bulk.clone_df, chunk_size = ideal_chunk_size) #(5)
	r1 = tabulate(clone_df1 = tr_search.clone_df, clone_df2 = tr_bulk.clone_df, pwmat = tr_search.rw_beta)
	outfile = os.path.join(dest, f"{sample_name}.{tag}.bulk_tabulation.tsv")
	print(f"WRITING: {outfile}")
	r1.to_csv(outfile, sep = '\t')
	toc = time.perf_counter()
	print(f"TABULATED IN {toc - tic:0.4f} seconds")
	del(tr_search)
	del(tr_bulk)
	return(f"{toc - tic:0.4f}s")


if __name__ == "__main__":
	# 888      888                               888             
	# 888      888                               888             
	# 888      888                               888             
	# 88888b.  888  8888b.     .d8888b   .d88b.  888888 .d8888b  
	# 888 "88b 888     "88b    88K      d8P  Y8b 888    88K      
	# 888  888 888 .d888888    "Y8888b. 88888888 888    "Y8888b. 
	# 888  888 888 888  888         X88 Y8b.     Y88b.       X88 
	# 888  888 888 "Y888888     88888P'  "Y8888   "Y888  88888P' 
	# December 18, 2020 
	# Seattle, WA
	from collections import defaultdict 
	import os 
	import pandas as pd
	import re
	# <path> path where bulk files can be found 
	path_bulkfiles = INSERT_HERE # <--------------------------

	# <path> where files reside
	path = os.path.join(path_to_base,'tcrdist','data','covid19')#os.path.join('data-raw','2020-08-31-mira_tcr_by_epitope/')
	# <all_files> list of all files
	all_files = [f for f in os.listdir(path) if f.endswith('.tcrdist3.csv')]
	# <restriction> list of tuples from Supporting Table S5: https://docs.google.com/spreadsheets/d/1WAmze6lir-v11odO-nYYbCiYVh8WhQh_vy2d1UPPKb0/edit#gid=942295061
	restriction = \
	[('m_55_524_ALR','A*01'),
	('m_1_8260_HTT','A*01'),
	('m_45_689_SYF','A*01'),
	('m_10_2274_LSP','B*07'),
	('m_155_59_RAR','B*07'),
	('m_133_102_NQK','B*15'),
	('m_48_610_YLQ','A*02'),
	('m_111_146_AEI','A*11'),
	('m_53_532_NYL','A*24'),
	('m_90_216_GYQ','C*07'),
	('m_140_92_NSS','A*01'),
	('m_55_524_ALR','B*08'),
	('m_183_39_RIR','A*03'),
	('m_10_2274_LSP','C*07'),
	('m_99_191_QEC','B*40'),
	('m_155_59_RAR','C*07'),
	('m_185_39_ASQ','B*15'),
	('m_147_73_DLF','B*08'),
	('m_110_148_ELI','B*44'),
	('m_51_546_AYK','A*03'),
	('m_44_697_FPP','B*35'),
	('m_118_136_ALN','A*11'),
	('m_176_44_SST','A*11'),
	('m_30_1064_KAY','B*57'),
	('m_192_31_FQP','B*15'),
	('m_70_345_DTD','A*01')]
	# <restrictions_dict> convert list to dictionary
	restrictions_dict = defaultdict(list)
	for k,v in restriction:
		restrictions_dict[k].append(v)
	# Loop through all files to construct a Dataframe of only those with strongest evidence of HLA-restriction
	cache = list()
	for f in all_files:
		rgs = re.search(pattern ='(mira_epitope)_([0-9]{1,3})_([0-9]{1,6})_([A-Z]{1,3})', 
			string = f).groups()	
		rgs4 = re.search(pattern ='(mira_epitope)_([0-9]{1,3})_([0-9]{1,6})_([A-Z]{1,4})', 
			string = f).groups()
		key3 = f'm_{rgs[1]}_{rgs[2]}_{rgs[3]}'
		key4 = f'm_{rgs4[1]}_{rgs4[2]}_{rgs4[3]}'
		setkey = f"M_{rgs[1]}"
		include = key3 in restrictions_dict.keys()
		alleles = restrictions_dict.get(key3)
		cache.append((setkey, key3, key4, f, int(rgs[1]), int(rgs[2]), include, alleles))

	# <hla_df>
	hla_df = pd.DataFrame(cache, columns = ['set', 'key3','key4','filename','set_rank','clones','hla_restricted','alleles']).\
		sort_values(['hla_restricted','clones'], ascending = True).\
		query('hla_restricted == True').reset_index(drop = True)

	for ind, row in hla_df.iterrows():
		print(ind)
		print(row)
		tag_level = '1E6'
		# <tag> the results with this symbol
		tag = f"{row['set']}_{tag_level}"
		ranked_centers_fn = f"{row['filename']}.ranked_centers_bkgd_ctlr_{tag_level}.tsv"
		benchmark_fn      = f"{row['filename']}.ranked_centers_bkgd_ctlr_{tag_level}.tsv.benchmark_tabulation.tsv"
		# <dest> place where .tsv files shall be saved
		dest = f'/fh/fast/gilbert_p/fg_data/tcrdist/t3/{tag}'
		if not os.path.isdir(dest):
			os.mkdir(dest)

		# <path> path where bulk files can be found 
		path = path_bulkfiles
		# <covid_smaples> list of all samples to consider>
		covid_samples = ["BS-EQ-09-T1-replacement_TCRB", "BS-EQ-0022-T0-replacement_TCRB", "860011133_TCRB", "BS-GIGI_86-replacement_TCRB", "BS-EQ-12-T1-replacement_TCRB", "BS-GIGI_44-replacement_TCRB", "KH20-09697_TCRB", "BS-GN-0007-T0-replacement_TCRB", "860011216_TCRB", "BS-EQ-0002-T1-replacement_TCRB", "KH20-09988_TCRB", "BS-GIGI_21-replacement_TCRB", "INCOV018-AC-5_TCRB", "860011205_TCRB", "BS-GN-06-T0-replacement_TCRB", "KHBR20-00107_TCRB", "550043156_TCRB", "550041451_TCRB", "KH20-09700_TCRB", "BS-EQ-25-T0-replacement_TCRB", "860011320_TCRB", "BS-GN-0008-T0-replacement_TCRB", "BS-EQ-0002-T2-replacement_TCRB", "860011277_TCRB", "KH20-11306_TCRB", "BS-EQ-10-T1-replacement_TCRB", "BS-EQ-18-T1-replacement_TCRB", "KH20-11638_TCRB", "860011493_TCRB", "BS-GN-0010-T0-replacement_TCRB", "860011220_TCRB", "KH20-09681_TCRB", "550041430_TCRB", "INCOV033-AC-3_TCRB", "KHBR20-00093_TCRB", "860011221_TCRB", "BS-GIGI_28-replacement_TCRB", "BS-GIGI_46-replacement_TCRB", "KH20-09665_TCRB", "BS-EQ-0018-T0-replacement_TCRB", "860011287_TCRB", "BS-GIGI_42-replacement_TCRB", "INCOV004-BL-5_TCRB", "KH20-11809_TCRB", "BS-EQ-0001-T1-replacement_TCRB", "KH20-09671_TCRB", "KHBR20-00111_TCRB", "BS-EQ-0027-T1-replacement_TCRB", "BS-EQ-0024-T0-replacement_TCRB", "BS-GN-13-T0-replacement_TCRB", "BS-EQ-11-T2-replacement_TCRB", "BS-GIGI_32-replacement_TCRB", "BS-EQ-31-T0-replacement_TCRB", "BS-EQ-0027-T2-replacement_TCRB", "BS-GIGI_06-replacement_TCRB", "KH20-09679_TCRB", "KH20-09676_TCRB", "860011462_TCRB", "BS-GIGI_56-replacement_TCRB", "KH20-09722_TCRB", "KH20-09751_TCRB", "550040014_TCRB", "BS-GIGI_96-replacement_TCRB", "BS-GN-0012-T0-replacement_TCRB", "KHBR20-00085_TCRB", "KH20-09667_TCRB", "INCOV086-AC-3_TCRB", "KH20-11307_TCRB", "BS-GIGI_49-replacement_TCRB", "INCOV020-AC-5_TCRB", "BS-GIGI_85-replacement_TCRB", "KH20-11650_TCRB", "BS-GIGI_15-replacement_TCRB", "BS-EQ-43-T0-replacement_TCRB", "BS-GIGI_50-replacement_TCRB", "INCOV013-AC-5_TCRB", "BS-EQ-0011-T1-replacement_TCRB", "BS-GIGI_66-replacement_TCRB", "860011140_TCRB", "BS-EQ-12-T2-replacement_TCRB", "KHBR20-00089_TCRB", "860011235_TCRB", "860011123_TCRB", "KHBR20-00121_TCRB", "KHBR20-00180_TCRB", "BS-GIGI_18-replacement_TCRB", "BS-GIGI_64-replacement_TCRB", "BS-EQ-0024-T1-replacement_TCRB", "550042607_TCRB", "KH20-11683_TCRB", "KHBR20-00179_TCRB", "BS-GN-0002-T0-replacement_TCRB", "550043906_TCRB", "BS-EQ-0021-T0-replacement_TCRB", "BS-EQ-07-T1-replacement_TCRB", "860011306_TCRB", "550043943_TCRB", "BS-GIGI_20-replacement_TCRB", "BS-GIGI_52-replacement_TCRB", "KH20-09970_TCRB", "860011129_TCRB", "KHBR20-00098_TCRB", "860011124_TCRB", "860011110_TCRB", "BS-GIGI_01-replacement_TCRB", "KH20-09670_TCRB", "BS-GIGI_72-replacement_TCRB", "BS-EQ-11-T3-replacement_TCRB", "860011120_TCRB", "INCOV008-BL-5_TCRB", "BS-GN-0004-T0-replacement_TCRB", "INCOV028-AC-3_TCRB", "KH20-11295_TCRB", "860011214_TCRB", "BS-GIGI_13-replacement_TCRB", "860011204_TCRB", "BS-GIGI_87-replacement_TCRB", "860011215_TCRB", "860011202_TCRB", "INCOV005-BL-5_TCRB", "BS-GIGI_10-replacement_TCRB", "860011489_TCRB", "550043905_TCRB", "BS-GN-01-T0-replacement_TCRB", "550042578_TCRB", "KHBR20-00118_TCRB", "KH20-09655_TCRB", "BS-EQ-0016-T1-replacement_TCRB", "KHBR20-00106_TCRB", "BS-GIGI_62-replacement_TCRB", "BS-EQ-37-T0-replacement_TCRB", "860011244_TCRB", "KHBR20-00157_TCRB", "BS-GIGI_80-replacement_TCRB", "KHBR20-00182_TCRB", "KHBR20-00206_TCRB", "860011318_TCRB", "KH20-11645_TCRB", "KH20-09673_TCRB", "550040140_TCRB", "BS-GIGI_36-replacement_TCRB", "KHBR20-00171_TCRB", "BS-GIGI_81-replacement_TCRB", "860011225_TCRB", "550043911_TCRB", "BS-GIGI_03-replacement_TCRB", "KH20-09698_TCRB", "BS-EQ-29-T0-replacement_TCRB", "860011239_TCRB", "KH20-11303_TCRB", "550041293_TCRB", "INCOV028-BL-3_TCRB", "INCOV004-AC-5_TCRB", "110047437_TCRB", "860011201_TCRB", "BS-EQ-23-T1-replacement_TCRB", "KHBR20-00119_TCRB", "KHBR20-00149_TCRB", "BS-EQ-0017-T0-replacement_TCRB", "KH20-11810_TCRB", "BS-GIGI_45-replacement_TCRB", "550042567_TCRB", "860011109_TCRB", "KH20-09752_TCRB", "KH20-09675_TCRB", "KH20-09987_TCRB", "BS-GIGI_41-replacement_TCRB", "BS-EQ-0008-T1-replacement_TCRB", "550042420_TCRB", "KH20-09747_TCRB", "KHBR20-00144_TCRB", "KH20-11675_TCRB", "INCOV031-BL-3_TCRB", "KH20-11697_TCRB", "KHBR20-00141_TCRB", "KH20-09948_TCRB", "INCOV036-AC-3_TCRB", "550041421_TCRB", "BS-EQ-0026-T0-replacement_TCRB", "BS-GIGI_08-replacement_TCRB", "BS-EQ-29-T1-replacement_TCRB", "KHBR20-00160_TCRB", "860011322_TCRB", "KHBR20-00165_TCRB", "860011327_TCRB", "KH20-09674_TCRB", "KH20-09661_TCRB", "BS-GIGI_16-replacement_TCRB", "BS-EQ-43-T1_BS-GIGI-137-replacement_TCRB", "860011117_TCRB", "INCOV033-BL-3_TCRB", "KHBR20-00148_TCRB", "BS-GIGI_79-replacement_TCRB", "550042515_TCRB", "550042526_TCRB", "INCOV003-BL-5_TCRB", "BS-EQ-35-T0-replacement_TCRB", "550042550_TCRB", "BS-GN-0011-T0-replacement_TCRB", "BS-EQ-07-T2-replacement_TCRB", "KHBR20-00103_TCRB", "INCOV030-AC-3_TCRB", "KH20-09997_TCRB", "BS-GIGI_61-replacement_TCRB", "KHBR20-00161_TCRB", "KH20-09996_TCRB", "INCOV030-BL-3_TCRB", "KHBR20-00156_TCRB", "KHBR20-00137_TCRB", "550041125_TCRB", "KHBR20-00130_TCRB", "860011498_TCRB", "KH20-11698_TCRB", "BS-EQ-17-T1b-replacement_TCRB", "KHBR20-00105_TCRB", "550043986_TCRB", "INCOV071-BL-3_TCRB", "KH20-11304_TCRB", "BS-EQ-0023-T0-replacement_TCRB", "KH20-11692_TCRB", "550042377_TCRB", "KH20-09720_TCRB", "BS-GIGI_35-replacement_TCRB", "KHBR20-00113_TCRB", "860011105_TCRB", "KH20-11294_TCRB", "KHBR20-00167_TCRB", "INCOV031-AC-3_TCRB", "KH20-09666_TCRB", "KHBR20-00081_TCRB", "BS-GIGI_70-replacement_TCRB", "BS-GIGI_55-replacement_TCRB", "INCOV009-BL-5_TCRB", "BS-EQ-34-T0-replacement_TCRB", "BS-GIGI_34-replacement_TCRB", "KHBR20-00185_TCRB", "BS-EQ-33-T0-replacement_TCRB", "BS-GIGI_19-replacement_TCRB", "KHBR20-00126_TCRB", "550042351_TCRB", "KH20-11647_TCRB", "BS-GIGI_74-replacement_TCRB", "860011229_TCRB", "KH20-11640_TCRB", "BS-GN-0005-T0-replacement_TCRB", "KHBR20-00083_TCRB", "KH20-09735_TCRB", "550042361_TCRB", "KH20-09963_TCRB", "KHBR20-00166_TCRB", "KH20-11695_TCRB", "KHBR20-00096_TCRB", "BS-GIGI_07-replacement_TCRB", "KH20-11811_TCRB", "INCOV069-BL-3_TCRB", "BS-GIGI_05-replacement_TCRB", "BS-GIGI_40-replacement_TCRB", "INCOV053-AC-3_TCRB", "860011116_TCRB", "INCOV016-BL-5_TCRB", "860011112_TCRB", "INCOV034-BL-3_TCRB", "INCOV032-BL-3_TCRB", "KH20-09664_TCRB", "550041314_TCRB", "860011206_TCRB", "550042640_TCRB", "INCOV017-BL-5_TCRB", "BS-EQ-10-T2-replacement_TCRB", "KHBR20-00172_TCRB", "KH20-11813_TCRB", "550043940_TCRB", "550041390_TCRB", "INCOV015-BL-5_TCRB", "KH20-09684_TCRB", "BS-EQ-36-T0-replacement_TCRB", "KH20-11637_TCRB", "BS-EQ-39-T0-replacement_TCRB", "BS-EQ-16-T2-replacement_TCRB", "KHBR20-00153_TCRB", "BS-GIGI_92-replacement_TCRB", "BS-GIGI_91-replacement_TCRB", "860011492_TCRB", "KHBR20-00117_TCRB", "550041499_TCRB", "KHBR20-00181_TCRB", "KHBR20-00173_TCRB", "KHBR20-00168_TCRB", "KHBR20-00123_TCRB", "KH20-09945_TCRB", "BS-GN-0014-T0-replacement_TCRB", "550043912_TCRB", "KH20-09955_TCRB", "KH20-11298_TCRB", "550043995_TCRB", "KHBR20-00170_TCRB", "INCOV070-BL-3_TCRB", "KH20-11305_TCRB", "KHBR20-00203_TCRB", "KHBR20-00188_TCRB", "BS-GIGI_53-replacement_TCRB", "INCOV036-BL-3_TCRB", "KH20-09959_TCRB", "KHBR20-00145_TCRB", "BS-EQ-33-T1-replacement_TCRB", "KH20-09677_TCRB", "INCOV087-BL-3_TCRB", "INCOV084-BL-3_TCRB", "KH20-11642_TCRB", "KH20-09687_TCRB", "KH20-09724_TCRB", "KH20-11687_TCRB", "BS-EQ-14-T3-replacement_TCRB", "KH20-11688_TCRB", "INCOV016-AC-5_TCRB", "BS-EQ-38-T1_BS-GIGI-117-replacement_TCRB", "INCOV075-AC-3_TCRB", "KHBR20-00162_TCRB", "550043973_TCRB", "550042259_TCRB", "550042183_TCRB", "BS-EQ-18-T2-replacement_TCRB", "KH20-11672_TCRB", "KH20-11643_TCRB", "KH20-11639_TCRB", "BS-GIGI_26-replacement_TCRB", "KH20-09696_TCRB", "KH20-09952_TCRB", "860011325_TCRB", "550043955_TCRB", "KHBR20-00177_TCRB", "BS-EQ-15-T1-replacement_TCRB", "BS-EQ-34-T1-replacement_TCRB", "860011312_TCRB", "KH20-11671_TCRB", "KHBR20-00196_TCRB", "KH20-09695_TCRB", "860011490_TCRB", "INCOV005-AC-5_TCRB", "550043980_TCRB", "KH20-09734_TCRB", "KH20-09750_TCRB", "KHBR20-00094_TCRB", "860011243_TCRB", "KHBR20-00143_TCRB", "BS-EQ-28-T1-replacement_TCRB", "KH20-09947_TCRB", "860011309_TCRB", "860011450_TCRB", "INCOV049-AC-3_TCRB", "860011310_TCRB", "860011113_TCRB", "KH20-09964_TCRB", "KH20-09991_TCRB", "KH20-11292_TCRB", "860011131_TCRB", "KHBR20-00200_TCRB", "KHBR20-00078_TCRB", "KHBR20-00186_TCRB", "BS-GIGI_27-replacement_TCRB", "550041250_TCRB", "KH20-09962_TCRB", "860011127_TCRB", "KH20-09999_TCRB", "860011338_TCRB", "KHBR20-00076_TCRB", "860011137_TCRB", "KHBR20-00204_TCRB", "550042523_TCRB", "INCOV017-AC-5_TCRB", "550043908_TCRB", "KHBR20-00205_TCRB", "KHBR20-00158_TCRB", "KH20-09708_TCRB", "860011218_TCRB", "860011382_TCRB", "BS-EQ-36-T1-replacement_TCRB", "860011125_TCRB", "BS-EQ-41-T1-replacement_TCRB", "BS-EQ-33-T2-replacement_TCRB", "860011132_TCRB", "550043971_TCRB", "860011246_TCRB", "INCOV006-BL-5_TCRB", "KHBR20-00134_TCRB", "550041349_TCRB", "KH20-10000_TCRB", "KH20-09718_TCRB", "KH20-09694_TCRB", "550042609_TCRB", "860011383_TCRB", "KH20-09741_TCRB", "550043927_TCRB", "KH20-11807_TCRB", "KH20-09668_TCRB", "KH20-11682_TCRB", "550042599_TCRB", "KHBR20-00127_TCRB", "550042631_TCRB", "KHBR20-00102_TCRB", "KH20-11814_TCRB", "KH20-09992_TCRB", "KHBR20-00178_TCRB", "KH20-09701_TCRB", "860011240_TCRB", "BS-EQ-41-T0-replacement_TCRB", "INCOV014-BL-5_TCRB", "KH20-11659_TCRB", "INCOV007-BL-5_TCRB", "KH20-09685_TCRB", "BS-EQ-28-T2-replacement_TCRB", "INCOV035-AC-3_TCRB", "INCOV011-AC-5_TCRB", "BS-GN-0016-T0-replacement_TCRB", "KHBR20-00152_TCRB", "BS-GN-0009-T0-replacement_TCRB", "INCOV025-AC-3_TCRB", "KHBR20-00128_TCRB", "550043167_TCRB", "INCOV051-AC-3_TCRB", "860011314_TCRB", "KHBR20-00108_TCRB", "KH20-09704_TCRB", "KH20-09749_TCRB", "KH20-09956_TCRB", "550042447_TCRB", "860011212_TCRB", "KH20-11670_TCRB", "BS-EQ-0028-T0-replacement_TCRB", "550040030_TCRB", "KH20-09746_TCRB", "KH20-11296_TCRB", "860011233_TCRB", "860011139_TCRB", "BS-GN-0003-T0-replacement_TCRB", "INCOV073-BL-3_TCRB", "860011203_TCRB", "KH20-09971_TCRB", "INCOV007-AC-5_TCRB", "KHBR20-00135_TCRB", "KH20-11291_TCRB", "860011241_TCRB", "KH20-09736_TCRB", "KHBR20-00201_TCRB", "BS-EQ-26-T1-replacement_TCRB", "KH20-11657_TCRB", "860011231_TCRB", "860011227_TCRB", "KHBR20-00082_TCRB", "KH20-09657_TCRB", "KH20-09966_TCRB", "KH20-09725_TCRB", "KH20-11681_TCRB", "550042442_TCRB", "KHBR20-00109_TCRB", "860011210_TCRB", "KH20-09985_TCRB", "860011238_TCRB", "550044028_TCRB", "860011445_TCRB", "550042400_TCRB", "BS-EQ-30-T2-replacement_TCRB", "KHBR20-00084_TCRB", "550042395_TCRB", "KHBR20-00136_TCRB", "KH20-11641_TCRB", "KH20-11302_TCRB", "KH20-11676_TCRB", "BS-EQ-0014-T1-replacement_TCRB", "KH20-11673_TCRB", "860011111_TCRB", "KH20-09967_TCRB", "KH20-09984_TCRB", "KH20-11300_TCRB", "KH20-11297_TCRB", "860011228_TCRB", "KH20-11301_TCRB", "KH20-11646_TCRB", "KHBR20-00151_TCRB", "KH20-11661_TCRB", "KH20-11299_TCRB", "KHBR20-00131_TCRB", "KH20-09709_TCRB", "860011496_TCRB", "KH20-09946_TCRB", "KHBR20-00163_TCRB", "550041393_TCRB", "KH20-09690_TCRB", "KHBR20-00080_TCRB", "KH20-09753_TCRB", "KH20-09669_TCRB", "KHBR20-00139_TCRB", "KH20-09691_TCRB", "INCOV019-AC-5_TCRB", "860011217_TCRB", "KH20-09980_TCRB", "INCOV034-AC-3_TCRB", "KH20-09951_TCRB", "BS-EQ-0003-T1-replacement_TCRB", "KH20-09978_TCRB", "KH20-11660_TCRB", "KH20-09961_TCRB", "INCOV035-BL-3_TCRB", "KH20-09990_TCRB", "860011303_TCRB", "KHBR20-00184_TCRB", "KH20-09954_TCRB", "INCOV044-AC-3_TCRB", "KHBR20-00154_TCRB", "KHBR20-00183_TCRB", "KHBR20-00122_TCRB", "KH20-09993_TCRB", "KH20-09748_TCRB", "550042547_TCRB", "KH20-09689_TCRB", "KH20-09727_TCRB", "BS-EQ-0014-T2-replacement_TCRB", "860011323_TCRB", "KH20-11812_TCRB", "KHBR20-00091_TCRB", "INCOV048-BL-3_TCRB", "KH20-09977_TCRB", "KH20-09986_TCRB", "KH20-09745_TCRB", "860011115_TCRB", "KHBR20-00087_TCRB", "550043914_TCRB", "550043981_TCRB", "860011499_TCRB", "INCOV054-BL-3_TCRB", "KHBR20-00100_TCRB", "550042239_TCRB", "KH20-11649_TCRB", "KH20-09953_TCRB", "860011304_TCRB", "KHBR20-00133_TCRB", "550042210_TCRB", "INCOV051-BL-3_TCRB", "860011497_TCRB", "KHBR20-00116_TCRB", "KH20-09715_TCRB", "BS-EQ-30-T1-replacement_TCRB", "860011348_TCRB", "860011118_TCRB", "KH20-09975_TCRB", "550043928_TCRB", "KH20-09714_TCRB", "860011242_TCRB", "KHBR20-00079_TCRB", "550041441_TCRB", "KH20-09960_TCRB", "860011326_TCRB", "KHBR20-00191_TCRB", "BS-GIGI_75-replacement_TCRB", "KH20-09968_TCRB", "KH20-09723_TCRB", "KH20-09979_TCRB", "KH20-11658_TCRB", "INCOV047-BL-3_TCRB", "INCOV002-AC-5_TCRB", "KH20-09731_TCRB", "INCOV029-AC-3_TCRB", "KH20-11674_TCRB", "110047542_TCRB", "KH20-09743_TCRB", "KH20-11678_TCRB", "KH20-11689_TCRB", "KH20-09699_TCRB", "KH20-09726_TCRB", "860011336_TCRB", "860011256_TCRB", "KHBR20-00097_TCRB", "INCOV006-AC-5_TCRB", "KHBR20-00086_TCRB", "KH20-11651_TCRB", "KHBR20-00125_TCRB", "KH20-09755_TCRB", "KH20-11655_TCRB", "KH20-09683_TCRB", "KH20-09976_TCRB", "860011280_TCRB", "INCOV047-AC-3_TCRB", "KH20-11684_TCRB", "860011219_TCRB", "KH20-09654_TCRB", "INCOV044-BL-3_TCRB", "KH20-09738_TCRB", "INCOV015-AC-5_TCRB", "INCOV050-AC-3_TCRB", "860011494_TCRB", "KH20-11662_TCRB", "KH20-11665_TCRB", "INCOV062-BL-3_TCRB", "860011226_TCRB", "KH20-09702_TCRB", "860011354_TCRB", "860011236_TCRB", "KHBR20-00095_TCRB", "860011209_TCRB", "860011317_TCRB", "KH20-09711_TCRB", "KH20-11699_TCRB", "KH20-09754_TCRB", "KHBR20-00190_TCRB", "860011388_TCRB", "860011356_TCRB", "550043920_TCRB", "KHBR20-00099_TCRB", "KHBR20-00129_TCRB", "KH20-11696_TCRB", "KH20-11686_TCRB", "KH20-11664_TCRB", "INCOV026-BL-3_TCRB", "KH20-11667_TCRB", "KH20-11666_TCRB", "860011347_TCRB", "KH20-11648_TCRB", "KH20-09739_TCRB", "550044000_TCRB", "KH20-09972_TCRB", "KH20-11668_TCRB", "860011138_TCRB", "860011208_TCRB", "550041239-1_TCRB", "550042542_TCRB", "860011223_TCRB", "KHBR20-00092_TCRB", "INCOV052-BL-3_TCRB", "INCOV050-BL-3_TCRB", "860011311_TCRB", "550044029_TCRB", "INCOV003-AC-5_TCRB", "KHBR20-00088_TCRB", "860011264_TCRB", "INCOV012-AC-5_TCRB", "INCOV001-AC-5_TCRB", "INCOV002-BL-5_TCRB", "KH20-11654_TCRB", "KHBR20-00115_TCRB", "KH20-09707_TCRB", "860011392_TCRB", "KH20-11656_TCRB", "KH20-09717_TCRB", "860011284_TCRB", "KHBR20-00174_TCRB", "KH20-11653_TCRB", "KH20-09719_TCRB", "KH20-09958_TCRB", "INCOV026-AC-3_TCRB", "860011119_TCRB", "KH20-09981_TCRB", "KH20-11652_TCRB", "860011250_TCRB", "KH20-09995_TCRB", "INCOV059-AC-3_TCRB", "860011282_TCRB", "KH20-11669_TCRB", "550041173_TCRB", "KH20-09973_TCRB", "KH20-09706_TCRB", "860011122_TCRB", "KH20-09969_TCRB", "KH20-11293_TCRB", "INCOV062-AC-3_TCRB", "860011340_TCRB", "INCOV055-BL-3_TCRB", "860011224_TCRB", "KH20-09658_TCRB", "KHBR20-00104_TCRB", "860011211_TCRB", "KH20-11679_TCRB", "860011130_TCRB", "860011495_TCRB", "INCOV027-BL-3_TCRB", "KH20-09693_TCRB", "KHBR20-00146_TCRB", "860011106_TCRB", "KH20-09721_TCRB", "INCOV077-BL-3_TCRB", "KH20-09716_TCRB", "INCOV027-AC-3_TCRB", "860011121_TCRB", "KH20-11644_TCRB", "INCOV045-BL-3_TCRB", "860011330_TCRB", "INCOV023-AC-3_TCRB", "860011260_TCRB", "KH20-09994_TCRB", "KH20-09744_TCRB", "KH20-09712_TCRB", "860011488_TCRB", "550042451_TCRB", "INCOV078-BL-3_TCRB", "KHBR20-00164_TCRB"]
		# <max number of cpus to use>
		ncpus = 2
		# <df_search> dataframe containing metaclonotypes
		df_search = pd.read_csv(os.path.join('hla_restricted_meta_clonotypes', ranked_centers_fn), sep = "\t")
		if df_search.shape[0] > 0:
			df_search = df_search[['cdr3_b_aa','v_b_gene','j_b_gene','pgen','radius','regex']]
			# <all_files> list of all files in the project <path>
			all_files = [f for f in os.listdir(path) if f.endswith('v_max.tsv')]
			# add file size
			all_files = [ (f, f.replace('.tsv.tcrdist3.v_max.tsv',''), os.stat(os.path.join(path,f)).st_size) for f in all_files ]
			# sort ascending by file size, smaller files first
			all_files = sorted(all_files, key=lambda x: x[2])
			# <all_files_ms> is <all_files> subset to only the samples with 30 days of COVID19 diagnosis, and not including ADPATIVE-COVID19 training samples
			all_files_ms = [(x,y,z) for x,y,z in all_files if y in covid_samples] 
			ts_iterable = [f for f,sn,_  in all_files_ms]
			import parmap
			timings = parmap.map(do_search2, ts_iterable, df_search = df_search.copy(), dest = dest, tag = tag, pm_processes=6, pm_pbar = True)
			pd.DataFrame({'file':ts_iterable, 'seconds':timings}).to_csv(os.path.join('hla_restricted_meta_clonotypes', benchmark_fn), sep = "\t")


# 8888888888 888b    888 8888888b.  
# 888        8888b   888 888  "Y88b 
# 888        88888b  888 888    888 
# 8888888    888Y88b 888 888    888 
# 888        888 Y88b888 888    888 
# 888        888  Y88888 888    888 
# 888        888   Y8888 888  .d88P 
# 8888888888 888    Y888 8888888P"


from mofchecker import MOFChecker
import os
import pandas as pd


def check(path):
	results = []
	descriptors = ['name','formula','has_carbon','has_hydrogen','has_atomic_overlaps','has_overcoordinated_c','has_overcoordinated_n','has_overcoordinated_h','has_undercoordinated_c',
	'has_undercoordinated_n','has_undercoordinated_rare_earth','has_metal','has_lone_molecule','metal_number','positive_charge_from_linkers','negative_charge_from_linkers',
	'possible_charged_fused_ring','has_suspicious_terminal_oxo','has_undercoordinated_alkali_alkaline','has_geometrically_exposed_metal']

	output_file = 'check.csv'
	existing_files = []
	if os.path.exists(output_file):
		df = pd.read_csv(output_file, header=None)
		existing_files = set(df[0].astype(str).tolist())
	if not os.path.exists(output_file):
		df = pd.DataFrame(columns=descriptors)
		df.to_csv(output_file, index=False)
	for structure in path:
		str_name = structure.split('.cif')[0]
		if str_name not in existing_files:
			mofchecker = MOFChecker.from_cif(os.path.join(cif_path,structure), primitive=False)
			results = mofchecker.get_mof_descriptors(descriptors)
			df = pd.DataFrame([results])
			df.to_csv(output_file, mode='a', header=False, index=False)	

if __name__ == "__main__":
	cif_path = './cifs'
	cifs = os.listdir(cif_path)
	cifs.sort()
	check(cifs)   


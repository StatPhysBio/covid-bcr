from os import listdir
from os.path import isfile, join
import json
from utils import *

def filter_anns(annotations, header_dict):
    ann_dict = {'productive': [],
                'unproductive': []}
    for ann in annotations:
        if ann['invalid']:
            continue
        #  Skip if has shm indels.
        if ann['has_shm_indels'][0]:
            continue
        #  Skip sequences with N in the middle.
        if remove_N_buffer(ann['input_seqs'][0]).count("N") != 0:
            continue
        # Unproductive, only oof
        if not ann['in_frames'][0]:
            ann['unique_ids'][0] = header_dict[ann['unique_ids'][0]]
            ann_dict['unproductive'].append(ann)
        #  Productive, only in frame and no stops
        elif ann['in_frames'][0] and not ann['stops'][0]:
            ann['unique_ids'][0] = header_dict[ann['unique_ids'][0]]
            ann_dict['productive'].append(ann)
    return ann_dict

#  .yaml is json, but json isn't yaml
#  partis saves files which are actually
#  .json format but named .yaml
def filter_partis(in_yaml):
    print("in file",in_yaml)
    save_name = in_yaml.replace(".yaml","_filtered.json")
    patient = in_yaml.split("/")[-1].split("_")[0]
    fastas = get_files("/gscratch/stf/zachmon/covid/patient_sorted_fastas/",".fasta")
    fasta_file = [f for f in fastas if f.split("/")[-1].split(".")[0] == patient][0]
    _, headers = fasta_read(fasta_file)
    header_dict = {}
    for h in headers:
        header_dict[h.split("|")[0].replace(":","c")] =  h

    print("save file",save_name)
    if isfile(save_name):
        print("Already filtered")
        return
    print("going to filter")
    with open(in_yaml) as f:
        annotations = json.load(f)['events']
        filtered_anns = filter_annotations(annotations, header_dict)
    print("saving!")
    with open(save_name, 'w') as outfile:
        json.dump(filtered_anns, outfile)
    print("saved",save_name)

#  TODO Filter anns from abstar
#def pickle_json(in_json):
#    json_data = []
#    pickle_name = in_json.replace(".json",".pickle")
#    with open(in_json, 'r') as f:
#        for i, json_dict in enumerate(f):
#            if len(json_dict) < 10:
#                continue
#            if 'd_gene' not in json_dict:
#                continue
#            data_dict = json.loads(json_dict)
#            json_data.append(data_dict)
#    with open(pickle_name, 'wb') as handle:
#        pickle.dump(json_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
#    print("Finished " + in_json)

def main():
    import argparse
    parser = argparse.ArgumentParser(
            description='filter annotations: remove sequences with N (not in padding),\n'
                        'removed shm indel seqs, divide prod and unprod')
    parser.add_argument('--yaml', type=str, help='path to yaml file')
    parser.add_argument('--json', type=str, help='path to json file')
    args = parser.parse_args()
    if args.yaml:
        filter_partis(args.yaml)
    if args.json:
        pickle_json(args.json)

if __name__ == '__main__':
    main()

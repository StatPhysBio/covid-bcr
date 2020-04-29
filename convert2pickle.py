import yaml
import json
import pickle
from os import listdir
from os.path import isfile, join

def pickle_yaml(in_yaml):
    print(in_yaml)
    pickle_name = in_yaml.replace(".yaml",".pickle")
    if isfile(pickle_name):
        print("Already pickled")
        return
    with open(in_yaml) as file:
        # The FullLoader parameter handles the conversion from YAML
        # scalar values to Python the dictionary format
        partis_output = yaml.load(file, Loader=yaml.CLoader)
        data_we_care_about = partis_output["events"]
    with open(pickle_name, 'wb') as handle:
        pickle.dump(data_we_care_about, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print("Finished " + in_yaml)

def pickle_json(in_json):
    json_data = []
    pickle_name = in_json.replace(".json",".pickle")
    with open(in_json, 'r') as f:
        for i, json_dict in enumerate(f):
            if len(json_dict) < 10:
                continue
            if 'd_gene' not in json_dict:
                continue
            data_dict = json.loads(json_dict)
            json_data.append(data_dict)
    with open(pickle_name, 'wb') as handle:
        pickle.dump(json_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print("Finished " + in_json)

def main():
    import argparse
    parser = argparse.ArgumentParser(
            description='yaml to pickle converter.')
    parser.add_argument('--yaml', type=str, help='path to yaml file')
    parser.add_argument('--json', type=str, help='path to json file')
    args = parser.parse_args()
    if args.yaml:
        pickle_yaml(args.yaml)
    if args.json:
        pickle_json(args.json)

if __name__ == '__main__':
    main()

import h5py
import argparse
from collections import defaultdict, OrderedDict
import yaml
from tqdm import tqdm

helen_positional_features = defaultdict(list)
helen_positional_labels = defaultdict(int)
all_helen_positions = set()
all_medaka_positions = set()
_gap_ = '*'
_ref_gap_ = '#'
_read_sep_ = ' '
_alphabet_ = 'ACGT'
_extra_bases_ = 'NRYSWKMBDHV'
decoding = _gap_ + _alphabet_.lower() + _alphabet_.upper() + _read_sep_ + _extra_bases_
encoding = OrderedDict(((a, i) for i, a in enumerate(decoding)))
helen_decoding = {0: '*', 1: 'A', 2: 'C', 3: 'G', 4: 'T', 5: '*'}


def read_helen_h5py(hdf5_image, with_labels):
    hdf5_file = h5py.File(hdf5_image, 'r')

    image_dataset = hdf5_file['image']
    position_dataset = hdf5_file['position']
    index_dataset = hdf5_file['index']

    label_dataset = None
    if with_labels:
        label_dataset = hdf5_file['label']

    total_records = image_dataset.shape[0]
    print("PARSING HELEN DICTIONARY")
    for hdf5_index in tqdm(range(0, total_records), ncols=100):

        position = position_dataset[hdf5_index]
        index = index_dataset[hdf5_index]
        image = image_dataset[hdf5_index]

        if with_labels:
            label = label_dataset[hdf5_index]
            for l, p, i, f in zip(label, position, index, image):
                helen_positional_features[(p, i)] = f
                helen_positional_labels[(p, i)] = l
                all_helen_positions.add((p, i))
        else:
            for p, i, f in zip(position, index, image):
                if p < 0 or i < 0:
                    continue
                helen_positional_features[(p, i)] = f
                all_helen_positions.add((p, i))


def comapare_medaka_h5py(medaka_h5py, with_labels):
    hdf5_file = h5py.File(medaka_h5py, 'r')
    used_labels_list = list(yaml.load(hdf5_file['medaka_feature_decoding'][()]))

    if with_labels:
        used_labels_list = [(base, count) for used, base, count in used_labels_list if used == True]
        map_labels_list = list(yaml.load(hdf5_file['medaka_label_counts'][()]))

    samples = hdf5_file['samples']
    with tqdm(total=len(samples), leave=True, ncols=100) as progress_bar:
        for sample in samples:
            total_records = hdf5_file.get('/samples/'+sample+'/features').shape[0]

            features_dataset = hdf5_file.get('/samples/'+sample+'/features')
            positions_dataset = hdf5_file.get('/samples/'+sample+'/positions')

            labels_dataset = None

            if with_labels:
                labels_dataset = hdf5_file.get('/samples/'+sample+'/labels')

            for i in range(0, total_records):
                medaka_features = list(features_dataset[i])
                position = tuple(positions_dataset[i])

                # compare features
                if position not in helen_positional_features.keys():
                    # pass
                    print("POSITION NOT PRESENT IN HELEN: ", position)
                else:
                    # remove from helen's set
                    all_medaka_positions.add(position)
                    # compare encoded feature's values
                    helen_features = helen_positional_features[position]
                    equal = True
                    max_diff = max(abs(helen_features - medaka_features))
                    if max_diff > 0.0001:
                        equal = False
                    if not equal:
                        print("FEATURES DID NOT MATCH IN POSITION: ", position)
                        print("MEDAKA FEATURES: ", ',\t'.join(str(round(x, 3)) for x in medaka_features))
                        print("HELEN  FEATURES: ", ',\t'.join(str(round(x, 3)) for x in helen_features))

                    # compare the labels
                    if with_labels:
                        label, count = tuple(labels_dataset[i])
                        label_base = decoding[label]
                        medaka_label = label_base
                        helen_label = helen_decoding[helen_positional_labels[position]]
                        if helen_label != medaka_label:
                            print("LABELS DID NOT MATCH IN POSITION: ", position, "HELEN LABEL: ", helen_label,
                                  "MEDAKA LABEL", medaka_label)
            progress_bar.refresh()
            progress_bar.update(1)
        progress_bar.close()

    print("HELEN EXTRA POSITIONS: ", all_helen_positions.difference(all_medaka_positions))
    print("MEDAKA EXTRA POSITIONS: ", all_medaka_positions.difference(all_helen_positions))


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--helen_h5py",
        type=str,
        required=True,
        help="H5PY file generated by HELEN."
    )
    parser.add_argument(
        "--medaka_h5py",
        type=str,
        required=True,
        help="H5PY file generated by MEDAKA."
    )
    parser.add_argument(
        "--with_label",
        type=lambda x: (str(x).lower() == 'true' or str(x).lower() == '1'),
        default=False,
        help="If true then compare labels too."
    )

    FLAGS, unparsed = parser.parse_known_args()
    read_helen_h5py(FLAGS.helen_h5py, FLAGS.with_label)
    comapare_medaka_h5py(FLAGS.medaka_h5py, FLAGS.with_label)

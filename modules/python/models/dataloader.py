import os
import numpy as np
import pandas as pd
from torch.utils.data import Dataset
import torchvision.transforms as transforms
import h5py
import torch


class SequenceDataset(Dataset):
    """
    Arguments:
        A CSV file path
    """

    def __init__(self, csv_path, transform=None):
        data_frame = pd.read_csv(csv_path, header=None, dtype=str)
        assert data_frame[0].apply(lambda x: os.path.isfile(x.split(' ')[0])).all(), \
            "Some images referenced in the CSV file were not found"
        self.transform = transforms.Compose([transforms.ToTensor()])
        self.file_info = list(data_frame[0])

    def __getitem__(self, index):
        # load the image
        hdf5_image = self.file_info[index]
        label_decoder = {'A': 0, 'C': 1, 'G': 2, 'T': 3, '_': 4}
        hdf5_file = h5py.File(hdf5_image, 'r')
        image_dataset = hdf5_file['simpleWeight']
        label_dataset = hdf5_file['label']

        image = []
        for image_line in image_dataset:
            image.append(list(image_line))

        labels = []
        for label in label_dataset:
            labels.append(label_decoder[chr(label[0])])

        hdf5_file.close()

        image = torch.Tensor(image)
        # label = torch.(image).type(torch.DoubleStorage)
        label = np.array(labels, dtype=np.int)
        # label = torch.from_numpy(label).type(torch.LongStorage)

        return image, label

    def __len__(self):
        return len(self.file_info)

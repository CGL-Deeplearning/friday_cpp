import argparse
import h5py
import sys
sys.path.append(".")
import ImageAnalyzer as image_analyzer

def analyze_image(hdf5_file_path):
    hdf5_file = h5py.File(hdf5_file_path, 'r')
    image_dataset = hdf5_file['image']
    image = image_dataset
    print(image.shape)
    image_analyzer.analyze_np_array(image, image.shape[0], image.shape[1])


if __name__ == '__main__':
    '''
    Processes arguments and performs tasks to generate the pileup.
    '''
    parser = argparse.ArgumentParser()
    parser.register("type", "bool", lambda v: v.lower() == "true")
    parser.add_argument(
        "--hdf5_file",
        type=str,
        required=True,
        help="Bed file containing confident windows."
    )
    FLAGS, unparsed = parser.parse_known_args()
    analyze_image(FLAGS.hdf5_file)
import glob
import os
from PIL import Image
from numpy import zeros


class ElementaryUnitCellImage:
    def __init__(self, start_path="."):
        self.lista_plikow = list((glob.glob(os.path.join(start_path, "*.png"))))

    def normalize_pixel_value(self, input_file):
        file = Image.open(input_file)
        normalized_pixel = [abs((k / 255) - 1) for k in list(file.getdata(0))]
        image_size = file.size
        pixel_column = list(zeros(image_size[1]))
        for i in range(0, image_size[1]):
            pixel_column[i] = normalized_pixel[image_size[0] * i:(i + 1) * image_size[0]]
        return pixel_column

    def apply_for_every_image_in_directory(self):
        return [self.normalize_pixel_value(str(k)) for k in self.lista_plikow]

__author__ = 'szymag'

from PIL import Image

im = Image.open("1.jpg")
im.rotate(45).show()
ls = list(im.getdata(0))
print(ls)

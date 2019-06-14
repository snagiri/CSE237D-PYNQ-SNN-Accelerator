################################ README ########################################
# This file is used to initialize the network with trained weights.'image_names'
# consists of names of the images that are needed to be read.
################################################################################

import numpy as np
import imageio

def learned_weights():
	image_names = ["1", "2", "3", "4", "5", "6"]
	ans = []
	for image in image_names:
		temp = []
		img = imageio.imread("training_images/" + image + ".png") 
		print("printing weights")
		#if(image == "1"):
			#print(img.shape)
		for i in img:
			#print("hohooh")
			#print(i)
			for j in i:
				if(j==0):
					temp.append(-0.7)
				else:
					temp.append(1)
		ans.append(temp)


	print("ansss")
	print(len(ans))
	print(len(ans[0]))
	fname = "wt"+".txt"
	np.savetxt(fname,ans,delimiter=' ',fmt = '%f')
	return ans

if __name__ == '__main__':
	a = learned_weights()
	print(a)

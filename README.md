# Gated Fusion Network for Single Image Dehazing

Wenqi Ren, Lin Ma, Jiawei Zhang, Jinshan Pan, Xiaochun Cao, Wei Liu, Ming-Hsuan Yang

In this paper, we propose an efficient algorithm to directly restore a clear image from a hazy input. The proposed algorithm hinges on an end-to-end trainable neural
network that consists of an encoder and a decoder. The
encoder is exploited to capture the context of the derived
input images, while the decoder is employed to estimate the
contribution of each input to the final dehazed result using the learned representations attributed to the encoder.
The constructed network adopts a novel fusion-based strategy which derives three inputs from an original hazy image by applying White Balance (WB), Contrast Enhancing
(CE) and Gamma Correction (GC). We compute pixel-wise
confidence maps based on the appearance differences between these different inputs to blend the information of the
derived inputs and preserve the regions with pleasant visibility. The final dehazed image is yielded by gating the important features of the derived inputs. To train the network,
we introduce a multi-scale based approach so that the halo
artifacts can be avoided. Extensive experimental results on
both synthetic and real-world images demonstrate that the
proposed algorithm performs favorably against the stateof-the-art algorithms.


# Scripts and pre-trained models
The pre-trained model can be found in models/snapshot.      
Users can use demo_test to generate the results of any images.  


# Citations
Please cite this paper in your publications if it helps your research:    
@inproceedings{Ren-CVPR-2018,    
   author = {Ren, Wenqi and Ma, Lin and Zhang, Jiawei and Pan, Jinshan and Cao, Xiaochun and Liu, Wei and Yang, Ming-Hsuan},    
   title = {Gated Fusion network for Single Image Dehazing},    
   booktitle = {IEEE Conference on Computer Vision and Pattern Recognition},   
   year = {2018}
}

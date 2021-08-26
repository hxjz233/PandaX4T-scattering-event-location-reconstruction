# R-L + BGMM

## input
txt file under dir.
### format: 
"<PMTID> <x_coordinate> <y_coordinate> <intensity>" for each line
last 8 lines:
-8:-3 : 3 pairs of other algorithms' results
-2:-1 : z coordinate evaluated, energy evaluated
The 8 lines are not necessary for the algorithm itself, but for comparison and visualization
                                                
## output:
graph showing: 
1. direct interpolation of PMT intensities;
2. image refering to signal processed by R-L algorithm;
3. BGMM Clustering result of deconvoluted image;
4. Comparison of: other 3 algorithms, the last BGMM clustering result, and the BGMM clustering result of direct interpolation image, which took event amount estimated by last result as input.

## intermediate variables:
- PSF: Take single scattering events and evaluate their signal COG. For events at rCOG <= 300mm, stack up the signal, each COG aligned to center, relative intensity to max signal intensity as values. Add the signals up and chop it to size 581mm * 581mm, finally unify the whole spread function. In practice, we only use the central 401mm * 401mm as PSF.
    
## Algorithm outline:
    **Take input;**
    
    **Interpolate signal on 1200mm x 1200mm meshgrid.** For points outside farthest PMT(cannot be interpolated), fill with intensity 10(somehow may change according to different event energy scale. The results are not greatly affected if filled values are at the same scale with signal background).
    
    **Estimate approximate event location.**
        **Estimate largest interpolated intensity location.**
        **if rLargest >= 560mm**
            **event location is at 570mm**
            (In earlier tests, max intensity location shows at 570mm roughly for any event at 540mm - 600mm radial distance. Have not come up with a valid algorithm to evaluate these event yet.)
        else
            use COG to estimate approximate location.
            
    if rapproximate <= 300mm (central events)
        **R-L deconvolution**
    else if 300mm < rapproximate < 560mm (margin events)
        **rotate the signal image** so that the event aligns to \theta = 0. 
        **mirror the signal** by line x = 600mm.
        **linearly interpolate** signals between.(In earlier tests, cubic interpolation raises more noise)
        **R-L deconvolution**
        
    **Cluster the deconvoluted signal.**
        for every single intensity on the meshgrid, change the value to its relative intensity to the largest intensity on the gird. Then multiply them by a weight, serving as threshold of collecting the specific intensity or not, and by more times or less. For this code, we take the number 20, and collect the specific coordinate int(20 * relative intensity) times.
        (BGMM takes discrete point coordinates as input, and gives a way of clustering them binormally, while giving the centers a weight referring to the ratio of how many points belonged to it. It is necessary to turn float numbers referring to intensity on the mesh into discrete coordinate sets.
         Practically, no obvious change occured if a random number was added to the coordinates, dispersing the points.)
        if weight of any cluster set < multithreshold
            the component is abandoned. Event numbers should have been less than that set in forehand
        This steps gives a reference for exactly how many events occured.
        
    **Cluster the interpolated intensity**
        Similarly, point_weight changing float numbers into point numbers is necessary. It is set 2 for the case.
        (The earlier clustering, however, brings the center closer to each other than expected. It is supposed that instant cluster to the intensity signal may bring a more likely result.)
        For a central case, the cluster center number is set the same as the result of the last step, and clusters directly on the interpolated intensity.
        For a margin event, the cluster center number is set 2 * eventnum, and clusters on the whole mirrored and re-interpolated signal, considering that it may perform better if the event and 'mirrored event' signal were put together and clustered.
        
    **Visualization** of above variables (Rotate coordinates back if necessary)
    
## Possible arguments for optimization:
    - radial borderline deciding if a central/margin event happened, 300mm at the time
    - multithreshold, deciding when to abandon a cluster set
    - deconvolution iteration, which is proved to converge to the initial signal without any dispersion, but PSNR may rise with rising iteration at first, and decreces as iteration grows even larger. Have not gone through mathematical details for this yet.
    - pointweight, deciding the fineness as collecting intensity when need to cluster
    - cluster covariance, served as input for BGMM function. Preferred 'spherical' or 'full'(default)
    - cluster iteration, 2000 in this code is enough for clustering process to converge, but may be optimized if the function logged that it did not converge for certain input
    - PSF size, according to certain convolution process in R-L, smaller PSF may bring faster outputs. It had been tested that 201mm * 201mm PSF does not restore the signal well. It should be mentioned that the speed of the whole progress is decided by clustering progress and deconvolution progress together.
        
## Other issues:
    Deconvolution process may show 0/0 warnings and output NaN arrays. Negatives may show up at specific round of iteration, thus cause a 0 intermidiate value in further convolutions and zero divide exceptions(which belonged to detailed process in R-L algorithm). The reason however remained unknown, but the issue often show up in margin cases and around abrupt values.
    
    Some margin events may throw 3 or 1 event center results at the last cluster step, when 2 events were estimated forehand. This means the last cluster centers are not symmetric to x = 600mm on the whole mirrored signal, and these events should however be clustered again, or even with random initial cluster centers.

# PSF template calculation (first attempt)

## input
  Event data file with the same format in RL+BGMM.
  It is preferred to input a series of low-energy real run event data, so that the template would be mostly set up by single-scattering events, and adds up to the diffusion of a single event.
  
## output:
  "templatewhole.csv", at size 581 * 581, presenting average diffusion of an event
  
## Algorithm outline:
  **Take input**
  **Estimate the location of max intensity.**
    and only take Rmax <= 300mm events as filted input for PSF calculation.
  **Interpolate signal on 1200mm x 1200mm meshgrid.**
  **Use COG to estimate event center location.**
    It is possible to estimate event location with other algorithms, which may lead to different PSFs. However, we do suppose that COG estimation is enough for the central events.
  **Change the frame, so that the signals align to each other at their COGs**
  **Add up the signals and Unify the PSF, so that the sum of PSF meshgrid values adds to 1**
  **Save the PSF result to .csv format**

## Possible arguments for optimization:
  - PSFsize, which decides the final output PSF's size by side length = 2 * PSFsize + 1
  - Center alignment. In this case we used COG to estimate event center, but the validity of the method is not guaranteed, and may lead to system errors.

## Other issues:
      
# PSF template calculation (refinement)

## input & output
  The same as PSF template calculation (first attempt)
      
## Algorithm outline:
  The same as PSF template calculation (first attempt), despite that: for the filter, also used PSF derived from "first attempt" and RL+BGMM to estimate whether the event do be a single scattering event.

## Possible arguments for optimization:
  This part pretty much relies on the precision of PSF derived from "first attempt" & the validity of the RL+BGMM method, so there is no explicit arguments for the optimization.
      
## Other issues:

# Pure BGMM

@author: Zhuowei Zhang

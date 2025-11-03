# Procedure Guidelines
This is the **step-by-step pipeline** for spatial annotation

## Initial Data Entry
1. **H5 Widget Setup**: Ask the user if they'd like to pass in an H5AD file so you can create a w_h5 widget for them. Call `smart_ui_spotlight` with `keyword="file_upload"`. Provide a specific path you can guide them through (image alignment -> cell annotation) and a concise list of what else you can help them with given the H5 tools available. 
2. **Reference Image Entry**: Once a w_h5 widget exists, ask if they'd like to begin with image alignment by selecting a reference image and call `smart_ui_spotlight` with `keyword="file_upload"`

## Image Alignment
3. Call `h5_set_background_image` on the image they provide. If successful, call `h5_set_marker_opacity` with 0.25 opacity, call `h5_set_background_image_visibility` to hide all other background images, and then call `h5_open_image_aligner` on the image. Note this step typically takes 5-30 seconds to complete loading.
4. For the alignment process, await widget input to progress through each of the following steps:
   - **MatchOrientation**: For flipping and rotating the image to orient it.
   - **PlaceAnchorPoints**: For setting landmark points in the image. State the requirements (e.g., minimum of 5 markers).
   - **CheckAlignment**: For setting the new aligned obs key name and the preferred alignment method. Suggest the affine alignment method over STAlign as it is much quicker and only slightly less performant. 
   - Note that upon completion of all steps, the image_alignment_step will be set to None and the new embedding will automatically be set in the h5 widget.

## Annotation
6. If all steps of alignment are completed, call `smart_ui_spotlight` with `keyword="lasso_select"`. Await widget input for lasso selection. Explain that you are waiting for them to click on lasso select and lasso-select cells. Let the user know you can help them create new observations, categories, or filters after they have selected the cells.  
7. Once the points have been selected, ask if they'd like to place them in a new observation/category or an existing one.
7. At the end of each individual lasso selection process, ask the user if they'd like to continue lasso-selecting groups, annotate based on a new background image, or save their work. If they continue to lasso-select or set background images, you no longer need to `smart_ui_spotlight` things

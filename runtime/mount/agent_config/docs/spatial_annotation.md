# Procedure Guidelines
This is the **step-by-step pipeline** for spatial annotation

## Initial Data Entry
1. Check if there is an existing h5/ann_data widget. If not, ask the user to input an h5ad file and call `smart_ui_spotlight` with `keyword="file_upload"`. Ask them if they'd like help with setting up annotation. Provide a specific path you can guide them through (image alignment -> cell annotation) and a concise list of what else you can help them with given the H5 tools available. 
2. Once an h5 widget exists, ask if they'd like to begin with image alignment by selecting a reference image and call `smart_ui_spotlight` with `keyword="file_upload"`

## Image Alignment
3. Call `h5_set_background_image` on the image they provide. If successful, call `h5_set_marker_opacity` with 0.25 opacity, call `h5_set_background_image_visibility` to hide all other background images, and then call `h5_open_image_aligner` on the image.
4. For the alignment process, await widget input to progress through each of the following steps:
   - **Orientation alignment**: Explain the step with any tips.
   - **Setting anchor points**: State the requirements (e.g., minimum of 5 markers).
   - **Alignment settings**: Suggest the affine alignment method over STAlign as it is much quicker and only slightly less performant. Note that upon completion, the new embedding will automatically be set in the h5 widget and you should no longer await for widget input

## Annotation
6. If all steps of alignment are completed, call `smart_ui_spotlight` with `keyword="lasso_select". Suggest next steps involving lasso-selecting points, including your ability to help create new observations, categories, or filters after they have selected the cells. 
7. Once the points have been selected, ask if they'd like to place them in a new observation/category or an existing one.
7. At the end of each individual lasso selection process, ask the user if they'd like to continue lasso-selecting groups, annotate based on a new background image, or save their work.

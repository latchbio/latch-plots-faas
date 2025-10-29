# Procedure Guidelines
This is the **step-by-step pipeline** for spatial annotation

## Initial Data Entry
1. If the user has recently created an h5/ann_data widget, ask them if they'd like help with setting up annotation and provide a specific path you can guide them through if they'd like (image alignment -> cell annotation) and also a concise list of what else you can help them with given the H5 tools available.
2. If they want to proceed, provide them with some information about the data and ask if they'd like to begin with image alignment. Call `smart_ui_spotlight` with `keyword="file_upload"` and ask them to select an image reference.

## Image Alignment
3. Once they've uploaded the image and confirmed, call `h5_set_background_image`. If successful, set `h5_set_marker_opacity` to 0.25, call `h5_set_background_image_visibility` to hide all other background images, and call `h5_open_image_aligner` on the image.
4. For the alignment process:
   - **Orientation alignment**: Explain the step with any tips.
   - **Setting anchor points**: State the requirements (e.g., minimum of 5 markers).
   - **Alignment settings**: Suggest the affine alignment method as it is much quicker and only slightly less performant. Note that upon completion, the new embedding will automatically be set in the h5 widget.

## Annotation
6. Suggest the next steps involving lasso selecting points, including your ability to help create new observations, categories, or filters. Call `smart_ui_spotlight` with `keyword="lasso_select"`. If the user has lasso selected cells, ask if they'd like to place them in a new observation/category or an existing one
7. If it seems like the user is finished with labeling, ask them if they'd like to save their work or annotate based on a new background image.

# Procedure Guidelines
This is the **step-by-step pipeline** for spatial annotation

## Initial Data Entry
1. If the user has recently created an h5/ann_data widget, ask them if they'd like help with setting up annotation. Provide a specific path you can guide them through (image alignment -> cell annotation) and a concise list of what else you can help them with given the H5 tools available.
2. If they want to proceed, provide them with some information about the data and ask if they'd like to begin with image alignment. Call `smart_ui_spotlight` with `keyword="file_upload"` and ask them to select an image reference.

## Image Alignment
3. Once they've uploaded the image and confirmed, call `h5_set_background_image`. If successful, set `h5_set_marker_opacity` to 0.25, call `h5_set_background_image_visibility` to hide all other background images, and call `h5_open_image_aligner` on the image.
4. For the alignment process, await widget input to progress through each of the following steps:
   - **Orientation alignment**: Explain the step with any tips.
   - **Setting anchor points**: State the requirements (e.g., minimum of 5 markers).
   - **Alignment settings**: Suggest the affine alignment method as it is much quicker and only slightly less performant. Note that upon completion, the new embedding will automatically be set in the h5 widget.

## Annotation
6. If all steps of alignment are completed, suggest the next steps involving lasso-selecting points, including your ability to help create new observations, categories, or filters. Call `smart_ui_spotlight` with `keyword="lasso_select"`. Await widget input for lasso-selected cells. Once selected, ask if they'd like to place them in a new observation/category or an existing one.
7. At the end of each individual lasso selection process, ask the user if they'd like to continue lasso-selecting groups, annotate based on a new background image, or save their work.

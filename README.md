# Image Quilting
Image quilting takes a texture and samples a portion of it and quilts samples together to make a larger, seamless texture. This code contains three methods: random, border matching, and border matching with cutting.


## Random Quilting

The easiest option which samples random squares and puts them together to make a larger picture.

## Border Matching

The intermediate option which samples squares which share border properties calculated by the sum of squared differences done by template matching.

## Border Matching with Cutting

A step ahead of the intermediate option which samples squares which share border properties calculated by the sum of squared differences done by template matching. It then cuts off the parts of the sample which is dissimilar to its neighbors.

[![Image Quilting](https://res.cloudinary.com/marcomontalbano/image/upload/v1652983742/video_to_markdown/images/youtube--xPWD77vmT2U-c05b58ac6eb4c4700831b2b3070cd403.jpg)](https://www.youtube.com/watch?v=xPWD77vmT2U "Image Quilting")

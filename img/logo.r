# pak::pkg_install("hexSticker")

library(hexSticker)
sysfonts::font_add_google("Gochi Hand", "gochi")
## Automatically use showtext to render text for future devices

img <- magick::image_read("img/fig.svg")
logo <- magick::image_ggplot(img, interpolate = TRUE)

sticker(
  logo,
  package = "muchart",
  p_size = 18,
  s_width = 1,
  s_height = 1,
  s_x = 1,
  s_y = 0.84,
  h_fill = "#ffffff",
  h_color = "#ffcd05",
  p_color = "#000000",
  p_family = "gochi",
  h_size = 2.4,
  white_around_sticker = T,
  filename = "logo.png",
  url = "https://prdm0.github.io/muchart/",
  u_size = 4.5,
  spotlight = T,
  l_alpha = 0.6,
  dpi = 300,
  u_color = "#0F2536"
)

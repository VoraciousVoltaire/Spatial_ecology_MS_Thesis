# Representative colony images 

# Directory containing images
img_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/For_images/"

# List of image filenames
img_files <- list.files(img_dir, full.names = TRUE)

# Read images
img_list <- lapply(img_files, readPNG)

# Titles of the images
img_titles <- c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Iceland", "Faroe Islands", 
                "Jarsteinen", "Eynhallow", "Isle of Canna", "Inishkea", "Little Saltee")

# Plotting
n <- length(img_list)
n_rows <- 2
n_cols <- 5

# Create a new plot
par(mfrow = c(n_rows, n_cols))

# Plot images and add titles
for (i in 1:n) {
  plot(as.raster(img_list[[i]]), main = img_titles[i], col = terrain.colors(256))
}

# Reset par to default
par(mfrow = c(1, 1))

# Adding titles to this----

# Directory containing images
img_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/For_images/"

# List of image filenames
img_files <- list.files(img_dir, full.names = TRUE)

# Read images
img_list <- lapply(img_files, readPNG)

# Titles of the images
img_titles <- c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Iceland", "Faroe Islands", 
                "Jarsteinen", "Eynhallow", "Isle of Canna", "Inishkea", "Little Saltee")

# Plotting
n <- length(img_list)
n_rows <- 2
n_cols <- 5

# Create a new plot
par(mfrow = c(n_rows, n_cols), oma = c(2, 2, 2, 2))

# Plot images and add titles
for (i in 1:n) {
  plot(as.raster(img_list[[i]]), main = "", xlab = "", ylab = "", axes = FALSE)  # Plot image
  title(img_titles[i], line = -2)  # Add title on top
}

# Reset par to default
par(mfrow = c(1, 1), oma = c(0, 0, 0, 0))

# Increasing the size of each subplot now----

# Directory containing images
img_dir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/For_images/"

# List of image filenames
img_files <- list.files(img_dir, full.names = TRUE)

# Read images
img_list <- lapply(img_files, readPNG)

# Titles of the images
img_titles <- c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Iceland", "Faroe Islands", 
                "Jarsteinen", "Eynhallow", "Isle of Canna", "Inishkea", "Little Saltee")

# Plotting
n <- length(img_list)
n_rows <- 2
n_cols <- 5

# Increase the size of each subplot
par(mfrow = c(n_rows, n_cols), mai = c(0.05, 0.05, 0.05, 0.05))

# Plot images and add titles
for (i in 1:n) {
  plot(as.raster(img_list[[i]]), main = "", xlab = "", ylab = "", axes = FALSE)  # Plot image
  title(img_titles[i], line = -2, cex.main = 2)  # Add title on top
}

# Reset par to default
par(mfrow = c(1, 1), mai = c(1, 1, 1, 1))


output.path = paste(getwd(), "/output/", sep="")

if (!file.exists(path = output.path)) {
  dir.create(path = output.path, showWarnings = TRUE)
}

library("randomForestSRC")

# Veteran's Administration Lung Cancer Trial. Randomized trial of
# two treatment regimens for lung cancer.
data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time, status) ~ ., data = veteran, ntree = 100, importance="none")
pdf(paste(output.path, "VeteranError.pdf", sep=""), width=10, height=10)
plot(v.obj)
dev.off()

pdf(paste(output.path, "VeteranSurvival.pdf", sep=""), width=10, height=10)
plot.survival(v.obj)
dev.off()

# Women's Interagency HIV Study (WIHS).  Competing risk data set
# involving AIDS in women.
data(wihs, package = "randomForestSRC")
wihs.obj <- rfsrc(Surv(time, status) ~ ., wihs, nsplit = 3, ntree =100)
pdf(paste(output.path, "wihs.pdf", sep=""), width=10, height=10)
plot.competing.risk(wihs.obj)
dev.off()


# New York air quality measurements. Mean ozone in parts per billion.
airq.obj <- rfsrc(Ozone ~ ., data = airquality, na.action = "na.impute", importance="none")
pdf(paste(output.path, "airq.pdf", sep=""), width=10, height=10)
plot.variable(airq.obj)
dev.off()

# Edgar Anderson's iris data with morphologic variation of three
# related species.
iris.obj <- rfsrc(Species ~., data = iris, ntree = 100)
pdf(paste(output.path, "iris.pdf", sep=""), width=10, height=10)
## Plot the results.
plot(iris.obj)
dev.off()







sink(paste(output.path, "examples.txt", sep=""), append=FALSE)

# Survival Analysis with first and second order depths for all variables.
data(veteran, package = "randomForestSRC")
v.obj <- rfsrc(Surv(time, status) ~ . , data = veteran)
v.max <- max.subtree(v.obj)

# First and second order depths.
print(round(v.max$order, 3))

# The minimal depth is the first order depth.
print(round(v.max$order[, 1], 3))

# Strong variables have minimal depth less than or equal
# to the following threshold.
print(v.max$threshold)

# This corresponds to the set of variables.
print(v.max$topvars)

sink()



data(mtcars)
mtcars.mod <- mtcars
mtcars.mod$carb <- factor(mtcars.mod$carb)
mtcars.mod$cyl <- factor(mtcars.mod$cyl)
mtcars.mix <- rfsrc(Multivar(carb, mpg, cyl) ~ ., data = mtcars.mod)

plot.variable(mtcars.mix, outcome.target = 3, which.outcome = 1, partial = TRUE, nvar = 1)
plot.variable(mtcars.mix, outcome.target = 3, which.outcome = 2, partial = TRUE, nvar = 1)

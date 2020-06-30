context("Absolute risk plotting")

# Handling warning messages coming from montecarlo integration
handler_validmc <- function(msg) {
    if (any(grepl("out of range", msg))) invokeRestart("muffleWarning")
}

testthat::skip_if_not_installed("glmnet")


data("brcancer")
mod_cb_glm <- fitSmoothHazard(cens ~ estrec * log(time) +
                                  horTh +
                                  age +
                                  menostat +
                                  tsize +
                                  tgrade +
                                  pnodes +
                                  progrec,
                              data = brcancer,
                              time = "time", ratio = 1)

mod_cb_glmnet <- fitSmoothHazard(cens ~ estrec * log(time) +
                                  horTh +
                                  age +
                                  menostat +
                                  tsize +
                                  tgrade +
                                  pnodes +
                                  progrec,
                              data = brcancer,
                              time = "time", ratio = 1, family = "glmnet")

mod_cb_gam <- fitSmoothHazard(cens ~ estrec * log(time) +
                                  horTh +
                                  age +
                                  menostat +
                                  tsize +
                                  tgrade +
                                  pnodes +
                                  progrec,
                              data = brcancer,
                              time = "time", ratio = 1, family = "gam")

smooth_risk_glm <- absoluteRisk(object = mod_cb_glm,
                                newdata = brcancer[1:10, ])

smooth_risk_glmnet <- absoluteRisk(object = mod_cb_glmnet,
                                newdata = brcancer[1:10, ])

smooth_risk_gam <- absoluteRisk(object = mod_cb_gam,
                                newdata = brcancer[1:10, ])

test_that("no error in plot method for absRiskCB objects - ggplot", {

    outglm <- try(plot(smooth_risk_glm),
                  silent = TRUE)

    outglmnet <- try(plot(smooth_risk_glmnet),
                  silent = TRUE)

    outgam <- try(plot(smooth_risk_gam),
                  silent = TRUE)

    # specify id names
    outglm_names <- try(plot(smooth_risk_glm,
                             id.names = paste0("Covariate Profile ", 1:10),
                             legend.title = "Type",
                             xlab = "time (days)",
                             ylab = "Cumulative Incidence (%)"),
                        silent = TRUE)

    # not enough ID names supplied
    expect_warning(plot(smooth_risk_glm,
                        id.names = paste0("Covariate Profile ", 1:9),
                        legend.title = "Type",
                        xlab = "time (days)",
                        ylab = "Cumulative Incidence (%)"))

    expect_false(inherits(outglm, "try-error"))
    expect_false(inherits(outglmnet, "try-error"))
    expect_false(inherits(outgam, "try-error"))
    expect_false(inherits(outglm_names, "try-error"))
})


test_that("no error in plot method for absRiskCB objects - matplot", {

    outglm <- try(plot(smooth_risk_glm, gg = FALSE),
                  silent = TRUE)

    # adding to exisitng plot
    outglmnet <- try(plot(smooth_risk_glmnet, gg = FALSE, add = TRUE),
                     silent = TRUE)

    expect_false(inherits(outglm, "try-error"))
    expect_false(inherits(outglmnet, "try-error"))
})

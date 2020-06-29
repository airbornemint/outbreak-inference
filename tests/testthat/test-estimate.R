test_that("pspline.estimate.timeseries works", {
  import::from(mgcv, gam, s)
  import::from(magrittr, "%>%")
  import::from(dplyr, filter)
  import::from(plyr, join)

  cases = data.frame(
    time=seq(0, 51),
    cases=c(c(rep(1, 13), seq(1, 50, length.out=13), seq(50, 1, length.out=13), rep(1, 13)))
  )

  model = gam(cases ~ s(time, k=20, bs="cp", m=3), family=gaussian, data=cases)

  eps = 0.1
  estTimes = data.frame(time=seq(min(cases$time), max(cases$time), by=eps))

  estCases = pspline.estimate.timeseries(
    model, estTimes,
    pspline.outbreak.cases, samples=20, level=.95
  )

  testEstimates = cases %>% join(estCases, "time") %>% filter(cases > 5)
  testRange = range(testEstimates$cases.median / testEstimates$cases)
  expect_gt(testRange[1], 0.8)
  expect_lt(testRange[1], 1.25)
})

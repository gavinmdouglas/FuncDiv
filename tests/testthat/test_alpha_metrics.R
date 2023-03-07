# Test that all alpha diversity metrics are returning expected values.
# The expected values were primarily produced with skbio.diversity.alpha functions.
# The exceptions were for Faith's PD, which was checked against faith_pd in the abdiv package,
# the Gini-Simpson index, which was checked against the diversity function in the diverse R package,
# and the inverse Simpson index, which was checked against the diversity function in the vegan R package.

test_that("faiths_pd returns expected value", {
  expect_equal(faiths_pd(tips_in_sample = test_tree$tip.label[1:100],
                         tree = test_tree),
               48.471429)
})

test_that("richness returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "richness"), 6)
})

test_that("richness returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "richness"), 3)
})


test_that("shannon_index returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "shannon_index"), 2.374556, tolerance = 1e-7)
})

test_that("shannon_index returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "shannon_index"), 1.405639, tolerance = 1e-7)
})

test_that("berger_parker_dominance returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "berger_parker_dominance"), 0.294118, tolerance = 1e-5)
})

test_that("berger_parker_dominance returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "berger_parker_dominance"), 0.5)
})

test_that("ENS_pie returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "ENS_pie"), 4.737705, tolerance = 1e-7)
})

test_that("ENS_pie returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "ENS_pie"), 2.461538, tolerance = 1e-6)
})

test_that("fishers_alpha returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "fishers_alpha", min_unique = 3), 3.3050, tolerance = 1e-5)
})

test_that("fishers_alpha returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "fishers_alpha", min_unique = 3), 1.743444, tolerance = 1e-7)
})

test_that("heips_evenness returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "heips_evenness"), 0.837152, tolerance = 1e-6)
})

test_that("heips_evenness returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "heips_evenness"), 0.824676, tolerance = 1e-6)
})

test_that("margalefs_richness returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "margalefs_richness"), 1.7647806)
})

test_that("margalefs_richness returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "margalefs_richness"), 0.9617967)
})

test_that("mcintoshs_dominance returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "mcintoshs_dominance"), 0.713662)
})

test_that("mcintoshs_dominance returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "mcintoshs_dominance"), 0.56094742)
})

test_that("mcintoshs_evenness returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "mcintoshs_evenness"), 0.63984058)
})

test_that("mcintoshs_evenness returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "mcintoshs_evenness"), 0.8271702)
})

test_that("menhinicks_richness returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "menhinicks_richness"), 1.45521375)
})

test_that("menhinicks_richness returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "menhinicks_richness"), 1.06066017)
})

test_that("pielous_evenness returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "pielous_evenness"), 0.91860367)
})

test_that("pielous_evenness returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "pielous_evenness"), 0.8868595)
})

test_that("gini_simpson_index returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "gini_simpson_index"), 0.78892734)
})

test_that("gini_simpson_index returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "gini_simpson_index"), 0.59375)
})

test_that("simpsons_evenness returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "simpsons_evenness"), 0.78961749)
})

test_that("simpsons_evenness returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "simpsons_evenness"), 0.82051282)
})

test_that("inverse_simpson_index returns expected value", {
  expect_equal(compute_alpha_div(in_vec, "inverse_simpson_index"), 4.73770492)
})

test_that("inverse_simpson_index returns expected value, with missing data input", {
  expect_equal(compute_alpha_div(in_vec_complex, "inverse_simpson_index"), 2.46153846)
})

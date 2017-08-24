#!/usr/bin/env r
conf = as.real(argv[1])
data = read.table('input.data', header = TRUE)
result = t.test(data$x, data$y, conf.level=conf, alternative='greater')
# print(result)
# cat(sprintf("%f %f %f %f %f %f %f %f\n", result$statistic, result$parameter, result$p.value, result$conf.int[1], result$conf.int[2], attributes(result$conf.int)$conf.level, result$estimate[1], result$estimate[2]))
cat(sprintf("%e", result$p.value))

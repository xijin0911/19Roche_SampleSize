# 19Roche_SampleSize
This is the internship in Roche (Sep 2019- Nov 2019), about sample size calculation for prediction models.

## Main goal 

Xijin Chen worked in her internship to understand, evaluate, implement, compare, and present two methods to calculate the required sample size for multivariate prediction models. Prediction models are frequently requested by non-statistical stakeholders, e.g. to predict the occurrence of an adverse event based on patient characteristics. However, the sample size available from clinical trial datasets is frequently insufficient for this task. The goal of this internship was to provide guidance and tools to facilitate quick decision-making on whether or not the creation of a prediction model has any possibility of providing meaningful insights. The main references used were a theoretical approach in van Smeden et al (2018) and a simulation-based approach Riley et al (2018) for binary endpoints. Xijin’s main duties included:

•	Gained an understanding of the work being done within Roche Biostatistics, including the analysis of clinical trial data and the interpretation of clinical trial protocols. Familiarized herself with relevant Standard Operating Procedures and technical systems. Attended training events, statistical seminars, and Study Team Meetings as learning opportunities.
•	Reviewed the relevant literature and understood the approaches outlined in van Smeden et al (2018) and Riley et al (2018). Summarized the theoretical background, relevant definitions, and implementation of these approaches in a LaTeX document.

•	Investigated limitations of the approaches and the sensitivity of the calculated sample size on input parameters, such as Cox-Snell’s R². Provided guidance on how to correctly use these approaches and recommendations for input parameters.

•	Used simulation techniques in R to validate these approaches for typical scenarios at Roche and compared their results. Showed limitations e.g. when extrapolating beyond the data these methods were fitted on.

•	Created an R Shiny application that implements these methods to calculate sample size for a given scenario and provides guidance on their use to make it easy for members of the Biostatistics department to apply these methods in the future.

•	Summarized the internship results and presented them to the Biostatistics department, colleagues from Safety Science working on a related project, and at a meeting with another pharma company.

•	Created an R Shiny application that can be used by the Department to easily apply the methods of van Smeden et al (2018) and Riley et al (2018). This will be beneficial when engaging stakeholders on prediction models.

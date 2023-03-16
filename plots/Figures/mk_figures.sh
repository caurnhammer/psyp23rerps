# Check if tmp dir exists
[[ -d tmp ]] || mkdir tmp

# Copy predictor densities
convert -density 500 \
	../Stimuli/adbc23_Densities.pdf \
	-background white -alpha remove -alpha off \
	Stimuli_Densities.png

# Copy/create RT graphs
convert -density 500 \
	../rRTs_Plaus_Clozedist/RT_logRT.pdf \
	../rRTs_Plaus_Clozedist/RT_logRT_wavelegend.pdf \
	-append \
	-background white -alpha remove -alpha off \
	RT_logRT.png

convert -density 500 \
	../rRTs_Plaus_Clozedist/RT_est_PlausCloze_distractor.pdf \
	../rRTs_Plaus_Clozedist/RT_res_PlausCloze_distractor.pdf \
	+append -page 1800x900 \
	tmp/RT_est_res_PlausClozedist.pdf
convert	-density 2000 \
	tmp/RT_est_res_PlausClozedist.pdf \
	../rRTs_Plaus_Clozedist/RT_logRT_wavelegend.pdf \
	-append -page 400x250 \
	-background white -alpha remove -alpha off \
	RT_est_res_PlausClozedist.png

convert	-density 1000 \
	../rRTs_Plaus_Clozedist/RT_coefficients.pdf \
	../rRTs_Plaus_Clozedist/RT_zvalue.pdf \
	+append -page 400x200 \
	-background white -alpha remove -alpha off \
	RT_coefficients_zvalue.png

convert	-density 1000 \
	../rRTs_PrecritRT_Plaus_Clozedist/RT_coefficients.pdf \
	../rRTs_PrecritRT_Plaus_Clozedist/RT_zvalue.pdf \
	+append -page 400x200 \
	-background white -alpha remove -alpha off \
	RT_Precrit_coefficients_zvalue.png

# Copy/create ERP graphs
convert -density 500 \
	../rERPs_Plaus_Clozedist/Waveforms/Observed_Full.pdf \
	-background white -alpha remove -alpha off \
	ERP_Observed_Full.png

convert -density 500 \
	../rERPs_Plaus_Clozedist_across/Waveforms/t-values.pdf \
	-background white -alpha remove -alpha off \
	ERP_across_tvalues.png

convert -density 1000 \
	../rERPs_ReadingTime/Waveforms/Coefficients_Pz.pdf \
	-background white -alpha remove -alpha off \
	ERP_RT.png

convert -density 500 \
	../rERPs_Plaus_Clozedist/Waveforms/EstimatedPz_InterceptPlausCloze_distractor.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/ResidualPz_InterceptPlausCloze_distractor.pdf \
	+append -page 1680:840 \
	tmp/ERP_EstRes_Pz.pdf
convert -density 1000 \
	tmp/ERP_EstRes_Pz.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/EstimatedPz_InterceptPlausCloze_distractor_wavelegend.pdf \
	-append -page 800:500 \
	-background white -alpha remove -alpha off \
	ERP_EstRes_Pz.png

convert -density 500 \
	../rERPs_Plaus_Clozedist/Waveforms/Coefficients_C3.pdf  \
	../rERPs_Plaus_Clozedist/Waveforms/Coefficients_Pz.pdf \
	+append -page 1680:840 \
	tmp/ERP_Coef_C3Pz.pdf
convert -density 1000 \
	tmp/ERP_Coef_C3Pz.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/Coefficients_C3_wavelegend.pdf \
	-append -page 800:500 \
	-background white -alpha remove -alpha off \
	ERP_Coef_C3Pz.png

convert -density 500 \
	../rERPs_Plaus_Clozedist/Waveforms/EstimatedPz_InterceptCloze_distractor.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/EstimatedPz_InterceptPlaus.pdf \
	+append -page 1680:840 \
	tmp/ERP_Est_Iso_Pz.pdf
convert -density 1000 \
	tmp/ERP_Est_Iso_Pz.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/EstimatedPz_InterceptPlausCloze_distractor_wavelegend.pdf \
	-append -page 800:500 \
	-background white -alpha remove -alpha off \
	ERP_Est_Iso_Pz.png

magick -density 250 \
	../rERPs_Plaus_Clozedist/Topos/Observed_B_250-400.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_B_300-500.pdf \
    ../rERPS_Plaus_Clozedist/Topos/Observed_B_600-1000.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_topolegend.pdf \
	+append -page 400:130 \
	-background white -alpha remove -alpha off \
	ERP_Observed_Topos_B.png

convert -density 250 \
	../rERPs_Plaus_Clozedist/Topos/Observed_C_250-400.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_C_300-500.pdf \
    ../rERPs_Plaus_Clozedist/Topos/Observed_C_600-1000.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_topolegend.pdf \
	+append -page 400:130 \
	-background white -alpha remove -alpha off \
	ERP_Observed_Topos_C.png

convert -density 250 \
	../rERPs_Plaus_Clozedist/Topos/Estimated_InterceptCloze_distractor_B_600-1000.pdf \
    ../rERPs_Plaus_Clozedist/Topos/Estimated_InterceptPlaus_B_600-1000.pdf \
	../rERPs_Plaus_Clozedist/Topos/Estimated_InterceptPlausCloze_distractor_B_600-1000.pdf \
    ../rERPs_Plaus_Clozedist/Topos/Observed_B_600-1000.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_topolegend.pdf \
	+append -page 350:80 \
	-background white -alpha remove -alpha off \
	ERP_Estimated_Topos_B.png

convert -density 500 \
	multistreams.pdf \
	-background white -alpha remove -alpha off \
	multistreams.png

convert -density 500 \
	RI_model.pdf \
	-background white -alpha remove -alpha off \
	RI_model.png

# remove tmp
rm -r tmp
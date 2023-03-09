cp ../rERPs_Plaus_Clozedist/Waveforms/Observed_Full.pdf ERP_Observed_Full.pdf

montage -mode concatenate -density 500 \
	../rERPs_Plaus_Clozedist/Waveforms/EstimatedPz_InterceptPlausCloze_distractor.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/ResidualPz_InterceptPlausCloze_distractor.pdf \
	EstimatedPz_InterceptPlausCloze_distractor.pdf_wavelegend.pdf \
	-tile 2x2 \
	ERP_EstRes_Pz.pdf

montage -mode concatenate -density 500 \
	../rERPs_Plaus_Clozedist/Waveforms/Coefficients_C3.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/Coefficients_Pz.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/Coefficients_C3.pdf_wavelegend.pdf \
	-tile 2x2 \
	ERP_Coef_C3Pz.pdf

montage -mode concatenate -density 500 \
	../rERPs_Plaus_Clozedist/Waveforms/EstimatedPz_InterceptCloze_distractor.pdf \
	../rERPs_Plaus_Clozedist/Waveforms/EstimatedPz_InterceptPlaus.pdf \
	EstimatedPz_InterceptPlausCloze_distractor.pdf_wavelegend.pdf \
	-tile 2x2 \
	ERP_Est_Iso_Pz.pdf

convert -density 250 ../rERPs_Plaus_Clozedist/Topos/Observed_B_250-400.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_B_300-500.pdf \
    ../rERPS_Plaus_Clozedist/Topos/Observed_B_600-1000.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_topolegend.pdf \
	+append -page 400:130 \
	ERP_Observed_Topos_B.pdf

convert -density 250 ../rERPs_Plaus_Clozedist/Topos/Observed_C_250-400.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_C_300-500.pdf \
    ../rERPs_Plaus_Clozedist/Topos/Observed_C_600-1000.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_topolegend.pdf \
	+append -page 400:130 \
	ERP_Observed_Topos_C.pdf

convert -density 250 ../rERPs_Plaus_Clozedist/Topos/Estimated_InterceptCloze_distractor_B_600-1000.pdf \
    ../rERPs_Plaus_Clozedist/Topos/Estimated_InterceptPlaus_B_600-1000.pdf \
	../rERPs_Plaus_Clozedist/Topos/Estimated_InterceptPlausCloze_distractor_B_600-1000.pdf \
    ../rERPs_Plaus_Clozedist/Topos/Observed_B_600-1000.pdf \
	../rERPs_Plaus_Clozedist/Topos/Observed_topolegend.pdf \
	+append -page 350:80 \
	ERP_Estimated_Topos_B.pdf

cp ../rERPs_Plaus_Clozedist_across/Waveforms/t-values.pdf ERP_across_tvalues.pdf
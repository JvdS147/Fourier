make clean
export simulate_pattern_highest_peak="1000"
export simulate_pattern_background="30000"
export zero_point_error="0.06"
export preferred_orientation="0.7"
export full_width_half_max="0.25"
make
#cp Fourier Fourier-Simulate_Experimental_Pattern_${simulate_pattern_highest_peak}_${simulate_pattern_background}_${zero_point_error}_${preferred_orientation}_${full_width_half_max}

let INDEX=4
python3 calc_d2s/create_d2s_jobs.py --data_input_path /30days/s4430291/Genomes_for_AFphylogeny_red_40_${INDEX} --data_output_path /90days/s4430291/Genomes_for_AFphylogeny_red_40_${INDEX}_D2S --temp T --submit T --dry_run F --index=${INDEX}

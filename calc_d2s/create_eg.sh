let INDEX=3
python3 calc_d2s/create_d2s_jobs.py --data_input_path /30days/s4430291/Genomes_for_AFphylogeny_red_40_${INDEX} --data_output_path /90days/s4430291/Genomes_for_AFphylogeny_red_40_${INDEX}_D2S --temp F --submit F --dry_run F --index=${INDEX}
# python3 calc_d2s/create_d2s_jobs.py --data_input_path /30days/s4430291/Genomes_for_AFphylogeny --data_output_path /90days/s4430291/Genomes_for_AFphylogeny --temp F --submit F --dry_run F --index=${INDEX}

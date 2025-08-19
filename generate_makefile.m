function generate_makefile(num_cpu)
% write_makefile([num_cpu])
% Generate a Makefile in the current directory.
% num_cpu = number of parallel num_cpu for "make -j" (default = 6)

    if nargin < 1 || isempty(num_cpu)
        num_cpu = 6;
    end

    mf = fullfile('.', 'Makefile');
    fid = fopen(mf, 'w');
    if fid == -1
        error('Cannot open "%s" for writing.', mf);
    end

    % ---- File content ----
    fprintf(fid, '# Just type: make\n');
    fprintf(fid, '# (In an example directory or work directory)\n\n');
    fprintf(fid, '# This makefile does:\n');
    fprintf(fid, '# 1) Makes a sync of the distribution code to a Compile\n');
    fprintf(fid, '#    directory here (will create directory if not here yet)\n');
    fprintf(fid, '# 2) Copies local files in this directory into that Compile directory\n');
    fprintf(fid, '# 3) Compiles code into executable for the example.\n');
    fprintf(fid, '# 4) Copies ''roms'' from Compile dir into this dir.\n\n');

    fprintf(fid, '.SUFFIXES:\n\n');
    fprintf(fid, '.PHONY: all clean depend\n\n');

    fprintf(fid, 'all:\n');
    fprintf(fid, '#\tln -s $(ROMS_ROOT)/Examples/Neptune_input 2>/dev/null || :\n');
    fprintf(fid, '#\tln -s $(ROMS_ROOT)/Examples/input_data 2>/dev/null || :\n');
    fprintf(fid, '\trsync -a $(ROMS_ROOT)/src/*.F $(ROMS_ROOT)/src/*.h  Compile\n');
    fprintf(fid, '\tmake show_tag\n');
    fprintf(fid, '\t@rsync -a $(ROMS_ROOT)/src/*.opt Compile\n');
    fprintf(fid, '\t@rsync -a $(ROMS_ROOT)/src/Makedefs.inc Compile\n');
    fprintf(fid, '\t@rsync -a $(ROMS_ROOT)/src/Make.depend Compile\n');
    fprintf(fid, '\t@rsync -a $(ROMS_ROOT)/src/Makefile Compile\n');
    fprintf(fid, '\t@rm Compile/*.f 2>/dev/null || :\n');
    fprintf(fid, '\tcp -p *.h *.F *.opt Makedefs.inc Compile 2>/dev/null || :\n');
    fprintf(fid, '\tcd Compile; make depend 2>/dev/null || :\n');
    fprintf(fid, '\tcd Compile; make -j%d; mv roms ..\n\n', num_cpu);

    fprintf(fid, 'tag := $(shell git rev-parse HEAD)\n');
    fprintf(fid, 'tline1 := git_hash=\n');
    fprintf(fid, 'tline2:= ''git_hash=\\"$(tag)\\"''\n');
    fprintf(fid, 'show_tag:\n');
    fprintf(fid, '\t@sed -i ''/$(tline1)/c\\      $(tline2)'' Compile/add_git_hash.F\n\n');

    fprintf(fid, 'compile_clean:\n');
    fprintf(fid, '\trm -r Compile/ roms 2>/dev/null || :\n\n');

    fprintf(fid, 'work_clean:\n');
    fprintf(fid, '\trm -r Compile/ *.F *.h *.in *.sh *.nc roms 2>/dev/null || :\n\n');

    fprintf(fid, 'code_check_clean:\n');
    fprintf(fid, '\t@cd $(ROMS_ROOT)/Examples/Flux_frc        ; make compile_clean\n');
    fprintf(fid, '\t@cd $(ROMS_ROOT)/Examples/Pipes_ana       ; make compile_clean\n');
    fprintf(fid, '\t@cd $(ROMS_ROOT)/Examples/Pipes_real      ; make compile_clean\n');
    fprintf(fid, '\t@cd $(ROMS_ROOT)/Examples/Rivers_ana      ; make compile_clean\n');
    fprintf(fid, '\t@cd $(ROMS_ROOT)/Examples/Rivers_real     ; make compile_clean\n');
    fprintf(fid, '\t@cd $(ROMS_ROOT)/Examples/Tracers_passive ; make compile_clean\n');
    fprintf(fid, '\t@cd $(ROMS_ROOT)/Examples/WEC_real        ; make compile_clean\n');

    fclose(fid);
    fprintf('Makefile written in current directory with -j%d\n', num_cpu);
end

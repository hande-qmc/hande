OUT=fciqmc/hubbard_fciqmc.out fciqmc/hubbard_fciqmc_real.out ifciqmc/hubbard_ifciqmc.out ifciqmc/hubbard_ifciqmc_real.out

BLOCK=$(addsuffix .block,$(basename $(OUT))) fciqmc/hubbard_fciqmc_spawn_cutoff.block

REBLOCK=../../../../tools/reblock_hande.py 

all: $(BLOCK)

# Specific cases --- chose start point.
fciqmc/%.block: fciqmc/%.out
	$(REBLOCK) --quiet --start 30000 $< > $@

fciqmc/hubbard_fciqmc_spawn_cutoff.block: $(wildcard fciqmc/hubbard_fciqmc*.out)
	$(REBLOCK) --quiet --start 30000 $^ > $@


ifciqmc/%.block: ifciqmc/%.out
	$(REBLOCK) --quiet --start 30000 $< > $@

# Catch all --- block everything.
%.block: %.out
	$(REBLOCK) --quiet $< > $@

clean:
	rm -f $(BLOCK)

# Generate dual pairs of algebras corresponding to Galois
# representations in the LMFDB.

# Requires the following Python packages:
# - dual-pairs <https://gitlab.com/pbruin/dual-pairs>
# - lmfdb-lite <https://github.com/roed314/lmfdb-lite>

GP = gp -q -D parisizemax=4G
SAGE = sage
PYTHON = $(SAGE) --python3

# Evaluate $(1) in GP with arguments $(2), catching errors.
gp_eval = echo "print(iferr("$(1)",e,write(\"/dev/stderr\",e);quit(1)));" | $(GP) $(2)

BIG_IMAGE_LABELS = 8T12 8T23 24T201 24T1353

BIG_IMAGE_FILES = $(BIG_IMAGE_LABELS:=-fields.txt)

all: GL1-matches.txt $(patsubst %,GL2_F%-matches.txt, 2 3 5) GSp4_F2-matches.txt

.DELETE_ON_ERROR:

.PRECIOUS: %-fields.txt %-reps.txt %-traces.txt

.SECONDARY: %.dualpair %.traces

%-fields.txt:
	$(PYTHON) download_nf.py $* > $@

GL1-reps.txt:
	$(PYTHON) download_GL1_reps.py > $@

GL2_F%-reps.txt:
	$(PYTHON) download_reps.py $* > $@

GSp4_F2-reps.txt:
	$(PYTHON) download_reps.py 2 4 > $@

%.dualpair: | $(BIG_IMAGE_FILES)
	$(call gp_eval,"to_lmfdb_format(dual_pair_from_big_image_field(Polrev("$(shell sed -ne "s/\"$*\".\(.*\)/\1/p" $(BIG_IMAGE_FILES))")))",big_images.gp) > $@

%.traces: %.dualpair
	$(SAGE) -c 'from dual_pairs.dual_pair_import import dual_pair_import; print(dual_pair_import("'$<'").nice_model().frobenius_traces_lmfdb())' > $@

GL2_F%-traces.txt: $(BIG_IMAGE_FILES)
	$(MAKE) $(shell sed -e 's/^"\(.*\)".*/\1.traces/' $(BIG_IMAGE_FILES))
	(for f in $$(($* * $* - 1))*.traces; do echo '("'$$f'", '`cat $$f`')'; done) | sed -e s/\\.traces// > $@

GL1-matches.txt: GL1-reps.txt
	$(SAGE) -c 'load("GL1.sage"); match()' > $@

GL2_F2-matches.txt: GL2_F2-reps.txt
	$(SAGE) -c 'load("GL2_F2.sage"); match()' > $@

GL2_F%-matches.txt: GL2_F%-traces.txt GL2_F%-reps.txt
	$(SAGE) -c 'load("GL2_Fl.sage"); match('$*')' > $@

GSp4_F2-matches.txt: | GSp4_F2-reps.txt
	$(SAGE) -c 'load("GSp4_F2.sage"); match()' > $@

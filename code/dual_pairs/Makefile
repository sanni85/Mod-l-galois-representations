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

all: dim-1-matches.txt GL2_F3-matches.txt GL2_F5-matches.txt GSp4_F2-matches.txt

.DELETE_ON_ERROR:

.PRECIOUS: %-fields.txt %-reps.txt %-traces.txt

%-fields.txt:
	$(PYTHON) download_nf.py $* > $@

dim-1-reps.txt:
	$(PYTHON) download_dim1_reps.py > $@

GL2_F%-reps.txt:
	$(PYTHON) download_reps.py $* > $@

GSp4_F2-reps.txt:
	$(PYTHON) download_reps.py 2 4 > $@

big-image-pairs: $(BIG_IMAGE_FILES)
	$(MAKE) $(shell sed -ne 's/^"\(.*\)".*/\1.gp/p' $(BIG_IMAGE_FILES))

big-image-traces:
	$(MAKE) $(shell sed -e 's/^"\(.*\)".*/\1-traces.gp/' $(BIG_IMAGE_FILES))

%.gp: | $(BIG_IMAGE_FILES)
	$(call gp_eval,"to_lmfdb_format(dual_pair_from_big_image_field(Polrev("$(shell sed -ne "s/\"$*\".\(.*\)/\1/p" $(BIG_IMAGE_FILES))")))",big_images.gp) > $@

%-traces.gp: %.gp
	$(SAGE) -c 'from dual_pairs.dual_pair_import import dual_pair_import; print(dual_pair_import("'$<'").nice_model().frobenius_traces_lmfdb())' > $@

GL2_F%-traces.txt:
	$(MAKE) big-image-traces
	(for f in $$(($* * $* - 1))*-traces.gp; do echo '("'$$f'", '`cat $$f`')'; done) | sed -e s/-traces\\.gp// | sort -V > $@

dim-1-matches.txt: dim-1-reps.txt
	$(SAGE) -c 'load("dim1.sage"); match()' > $@

GL2_F%-matches.txt: GL2_F%-traces.txt GL2_F%-reps.txt
	$(SAGE) -c 'load("match.sage"); match('$*')' > $@

GSp4_F2-matches.txt: | GSp4_F2-reps.txt
	$(SAGE) -c 'load("GSp4_F2.sage"); match()' > $@

# Generate dual pairs of algebras corresponding to Galois
# representations with image contained in GL_2(F_3) or GL_2(F_5).

# Requires the following Python packages:
# - dual-pairs <https://gitlab.com/pbruin/dual-pairs>
# - lmfdb-lite <https://github.com/roed314/lmfdb-lite>

GP = gp -q -D parisizemax=4G
PYTHON = sage --python3

BIG_IMAGE_LABELS = 8T12 8T23 24T201 24T1353

BIG_IMAGE_FILES = $(BIG_IMAGE_LABELS:=-fields.txt)

all: dim-1-matches.txt GL2_F3-matches.txt GL2_F5-matches.txt GSp4_F2-matches.txt

.DELETE_ON_ERROR:

.PRECIOUS: %-fields.txt GL2_F%-reps.txt

%-fields.txt:
	$(PYTHON) download_nf.py $* > $@

GL2_F%-reps.txt:
	$(PYTHON) download_reps.py $* > $@

GSp4_F2-reps.txt:
	$(PYTHON) download_reps.py 2 4 > $@

GSp4_F2-matches.txt: | GSp4_F2-reps.txt
	sage -c 'load("GSp4_F2.sage")' > $@

dim-1-reps.txt:
	$(PYTHON) download_dim1_reps.py > $@

# Evaluate $(1) in GP with arguments $(2), redirecting standard output
# to $(3).
gp_eval_raw = echo $(1) | $(GP) $(2) > $(3)
gp_catch = "iferr("$(1)",e,write(\"/dev/stderr\",e);quit(1))"
# Evaluate $(1) in GP with arguments $(2), catching errors and writing
# the result to $(3).
gp_eval = $(call gp_eval_raw,"print("$(call gp_catch,$(1))");",$(2),$(3))

big-image-pairs: $(BIG_IMAGE_FILES)
	$(MAKE) $(shell sed -ne 's/^"\(.*\)".*/\1.gp/p' $(BIG_IMAGE_FILES))

big-image-traces:
	$(MAKE) $(shell sed -e 's/^"\(.*\)".*/\1-traces.gp/' $(BIG_IMAGE_FILES))

%.gp: | $(BIG_IMAGE_FILES)
	$(call gp_eval,"to_lmfdb_format(dual_pair_from_big_image_field(Polrev("$(shell sed -ne "s/\"$*\".\(.*\)/\1/p" $(BIG_IMAGE_FILES))")))",big_images.gp,$@)

GL2_F%-matches.txt: GL2_F%-traces.txt GL2_F%-reps.txt
	sage -c 'load("match.sage"); match('$*')' > $@

%-traces.gp: %.gp
	sage -c 'from dual_pairs.dual_pair_import import dual_pair_import; print(dual_pair_import("'$<'").nice_model().frobenius_traces_lmfdb())' > $@

GL2_F%-traces.txt:
	$(MAKE) big-image-traces
	(for f in $$(($* * $* - 1))*-traces.gp; do echo '("'$$f'", '`cat $$f`')'; done) | sed -e s/-traces\\.gp// | sort -V > $@

dim-1-matches.txt: dim-1-reps.txt
	sage -c 'load("dim1.sage"); match()' > $@

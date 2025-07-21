# Generate dual pairs of algebras corresponding to Galois
# representations with image contained in GL_2(F_3) or GL_2(F_5).

# Requires the following Python packages:
# - dual-pairs <https://gitlab.com/pbruin/dual-pairs>
# - lmfdb-lite <https://github.com/roed314/lmfdb-lite>

GP = gp -q -D parisizemax=4G
PYTHON = sage --python3

BIG_IMAGE_LABELS = 8T12 8T23 24T201 24T1353

BIG_IMAGE_FILES = $(BIG_IMAGE_LABELS:=-fields.txt)

all: matches.txt

.DELETE_ON_ERROR:

.PRECIOUS: %-fields.txt mod-%-reps.txt

%-fields.txt:
	$(PYTHON) download_nf.py $* > $@

mod-%-reps.txt:
	$(PYTHON) download_reps.py $* > $@

# Evaluate $(1) in GP with arguments $(2), redirecting standard output
# to $(3).
gp_eval_raw = echo $(1) | $(GP) $(2) > $(3)
gp_catch = "iferr("$(1)",e,write(\"/dev/stderr\",e);quit(1))"
# Evaluate $(1) in GP with arguments $(2), catching errors and writing
# the result to $(3).
gp_eval = $(call gp_eval_raw,"print("$(call gp_catch,$(1))");",$(2),$(3))

big-image-pairs: $(BIG_IMAGE_FILES)
	$(MAKE) $(shell sed -ne 's/^"\(.*\)".*/\1.gp/p' $(BIG_IMAGE_FILES))

%.gp: | $(BIG_IMAGE_FILES)
	$(call gp_eval,"to_lmfdb_format(dual_pair_from_big_image_field(Polrev("$(shell sed -ne "s/\"$*\".\(.*\)/\1/p" $(BIG_IMAGE_FILES))")))",big_images.gp,$@)

mod-%-matches.txt: mod-%-traces.txt mod-%-reps.txt
	sage -c 'load("match.sage"); match('$*')' | sort -V > $@

mod-%-traces.txt:
	$(MAKE) big-image-pairs
	sage -c 'load("traces.sage"); compute_traces('$*')' | sort -V > $@

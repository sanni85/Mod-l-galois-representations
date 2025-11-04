# Generate dual pairs of algebras corresponding to Galois
# representations in the LMFDB.

# Requires the following Python packages:
# - dual-pairs <https://gitlab.com/pbruin/dual-pairs>
# - lmfdb-lite <https://github.com/roed314/lmfdb-lite>

GP = gp -q -D parisizemax=4G
SAGE = sage
PYTHON = $(SAGE) --python3

# Evaluate $(1) in GP with arguments $(2), catching errors.
gp_eval = echo "iferr("$(1)",e,write(\"/dev/stderr\",e);quit(1))" | $(GP) $(2)

BIG_IMAGE_LABELS = 8T12 8T23 24T201 24T1353

BIG_IMAGE_FILES = $(BIG_IMAGE_LABELS:=-fields.tsv)

CODOMAINS = GL1 GL2_F2 GL2_F3 GL2_F5 GSp4_F2

all: dual-pairs.tsv

dual-pairs.tsv: $(CODOMAINS:=-dual-pairs.tsv)
	sort -V $^ > $@

check: $(CODOMAINS:%=check-%)

check-%: %-dual-pairs.tsv
	$(call gp_eval,"check_traces(\"$*\");",traces.gp)

.DELETE_ON_ERROR:

# Do not delete intermediate files
.SECONDARY:

%-fields.tsv:
	$(PYTHON) download_nf.py $* > $@

GL1-reps.tsv:
	$(PYTHON) download_GL1_reps.py > $@

GL2_F%-reps.tsv:
	$(PYTHON) download_reps.py $* > $@

GSp4_F2-reps.tsv:
	$(PYTHON) download_reps.py 2 4 > $@

%.dualpair: | $(BIG_IMAGE_FILES)
	$(call gp_eval,"to_lmfdb_format(dual_pair_from_big_image_field(Polrev("$(shell sed -ne "s/$*\t\(.*\)/\1/p" $(BIG_IMAGE_FILES))")))",big_images.gp) > $@

GL2_F%-traces.tsv: $(BIG_IMAGE_FILES)
	$(MAKE) `sed -ne "s/^\($$(($* * $* - 1)).*\)\t.*/\1.dualpair/p" $(BIG_IMAGE_FILES)`
	$(call gp_eval,"big_image_traces($*)",traces.gp) > $@

GL1-dual-pairs.tsv: GL1-reps.tsv
	$(call gp_eval,"match();",GL1.gp) > $@

GL2_F2-dual-pairs.tsv: GL2_F2-reps.tsv
	$(call gp_eval,"match();",GL2_F2.gp) > $@

GL2_F%-dual-pairs.tsv: GL2_F%-traces.tsv GL2_F%-reps.tsv
	$(SAGE) -c 'load("GL2_Fl.sage"); match('$*')' > $@

GSp4_F2-dual-pairs.tsv: GSp4_F2-reps.tsv
	$(call gp_eval,"match();",GSp4_F2.gp) > $@

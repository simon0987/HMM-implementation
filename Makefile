TARGETS=train test
MODEL=model_0?.txt
RESULT=result?.txt
RM=rm -f

all: $(TARGETS)

train: train_hmm.cpp hmm.h
	g++ $< -o $@
test: test_hmm.cpp hmm.h
	g++ $< -o $@

clean:
	$(RM) $(TARGETS) $(RESULT) $(MODEL)
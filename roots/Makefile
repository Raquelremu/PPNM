# Compiler settings
CSC = mcs
CFLAGS = -optimize+

# Targets
all: ExerciseA.exe

ExerciseA.exe: ExerciseA.cs
	$(CSC) $(CFLAGS) -out:$@ $<

run: ExerciseA.exe
	mono ExerciseA.exe

clean:
	rm -f ExerciseA.exe

.PHONY: all run clean

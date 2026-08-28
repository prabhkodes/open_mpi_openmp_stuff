TARGET = sphere_interp
all: $(TARGET)

include ../common.mk

$(TARGET): sphere_interp.c
	$(CC) $(CFLAGS) $< -o $@ $(LDLIBS)

.PHONY: run scaling polar
run: $(TARGET)
	./$(TARGET) 20000 20000

## scaling: fit the complexity exponent for each method (TODO 6a)
## Start small -- brute force at 640k is 4e11 evaluations.
scaling: $(TARGET)
	@for n in 5000 10000 20000 40000; do \
	  echo "---- n = $$n ----"; \
	  ./$(TARGET) $$n $$n | grep -E 'brute|bucket|k-d'; \
	done

## polar: same point count, but check the pole degradation (TODO 6b)
polar: $(TARGET)
	./$(TARGET) 50000 50000

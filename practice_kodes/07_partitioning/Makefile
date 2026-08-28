TARGET = partition
all: $(TARGET)

include ../common.mk

$(TARGET): partition.f90
	$(FC) $(FFLAGS) $< -o $@

.PHONY: run scaling picture
run: $(TARGET)
	./$(TARGET) 256 256 16

## scaling: edge cut vs partition count (TODO 5b)
scaling: $(TARGET)
	@for p in 4 16 64 256; do \
	  echo "==== nparts = $$p ===="; \
	  ./$(TARGET) 256 256 $$p | sed -n '5,12p'; \
	done

## picture: render the partition so you can eyeball the curve
picture: $(TARGET)
	./$(TARGET) 256 256 16 >/dev/null
	@echo "wrote partition.ppm"
	@command -v magick >/dev/null 2>&1 && magick partition.ppm partition.png \
	  && echo "converted to partition.png" \
	  || echo "(install imagemagick for PNG, or open the .ppm directly -- macOS Preview reads it)"

clean-img:
	rm -f partition.ppm partition.png

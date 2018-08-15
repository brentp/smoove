package lumpy

import "testing"

func TestStartEndFix(t *testing.T) {
	in := "1	1116265	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116258"
	out := fixStartEnd(in)
	if out != "1	1116258	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265" {
		t.Errorf("didn't switch")
	}

	in = "1	1116265	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116258;"
	if fixStartEnd(in) != "1	1116258	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265;" {
		t.Errorf("didn't switch")
	}

	in = "1	1116258	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265;"
	if fixStartEnd(in) != in {
		t.Errorf("unneeded switch")
	}

	in = "1	1116265	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=111;"
	if fixStartEnd(in) != "1	111	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265;" {
		t.Errorf("improper switch: %s", fixStartEnd(in))
	}

	in = "1	1116265	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=111"
	if fixStartEnd(in) != "1	111	3510	N	<DEL>	0	.	SVTYPE=DEL;SVLEN=7;END=1116265" {
		t.Errorf("improper switch: %s", fixStartEnd(in))
	}

}

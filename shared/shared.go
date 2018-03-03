package shared

import (
	"io"
	"log"
	"os"
	"os/exec"
	"regexp"
	"strings"

	"github.com/biogo/hts/bam"
	"github.com/pkg/errors"
)

type Logger struct {
	*log.Logger
}

const Prefix = "[lumpy-smoother]"

func (l *Logger) Write(b []byte) (int, error) {
	l.Logger.Printf(string(b))
	return len(b), nil
}

var Slogger *Logger

func init() {
	l := log.New(os.Stderr, "[lumpy-smoother] ", log.Ldate|log.Ltime)
	Slogger = &Logger{Logger: l}
}

func Contains(haystack []string, needle string) bool {
	for _, h := range haystack {
		if h[0] != '~' && h == needle {
			return true
		}
		if h[0] == '~' {
			if match, _ := regexp.MatchString(h[1:], needle); match {
				return true
			}
		}
	}
	return false
}

type reader struct {
	io.ReadCloser
	cmd *exec.Cmd
}

func (r *reader) Close() error {
	if err := r.cmd.Wait(); err != nil {
		return errors.Wrap(err, "error closing cram reader")
	}
	return r.ReadCloser.Close()
}

// NewReader returns a bam.Reader from any path that samtools can read.
func NewReader(path string, rd int, fasta string) (*bam.Reader, error) {
	var rdr io.Reader
	if strings.HasSuffix(path, ".bam") {
		var err error
		rdr, err = os.Open(path)
		if err != nil {
			return nil, err
		}

	} else {
		cmd := exec.Command("samtools", "view", "-T", fasta, "-u", path)
		cmd.Stderr = Slogger
		pipe, err := cmd.StdoutPipe()
		if err != nil {
			return nil, errors.Wrap(err, "error getting stdout for process")
		}
		if err = cmd.Start(); err != nil {
			pipe.Close()
			return nil, errors.Wrap(err, "error starting process")
		}
		rdr = &reader{ReadCloser: pipe, cmd: cmd}
	}
	return bam.NewReader(rdr, rd)
}

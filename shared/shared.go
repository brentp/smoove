package shared

import (
	"log"
	"os"
	"regexp"
)

type Logger struct {
	*log.Logger
}

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

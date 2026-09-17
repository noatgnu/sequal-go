package sequal

import (
	"fmt"
	"regexp"
)

// crosslinkLabelPattern matches a crosslink label (e.g. "XL1"), distinct from an ambiguity label.
var crosslinkLabelPattern = regexp.MustCompile(`^XL[A-Za-z0-9]+$`)

// labelKind classifies a #label by its text shape, not by its internal PipeValue type.
type labelKind int

const (
	labelKindAmbiguity labelKind = iota
	labelKindCrosslink
)

func classifyLabel(label string) labelKind {
	if crosslinkLabelPattern.MatchString(label) {
		return labelKindCrosslink
	}
	return labelKindAmbiguity
}

// validateAmbiguityLabels checks every ambiguity-group reference resolves to a definition.
func validateAmbiguityLabels(mods []*Modification) error {
	defined := make(map[string]bool)
	for _, m := range mods {
		label := m.GetCrosslinkID()
		if label == nil || classifyLabel(*label) != labelKindAmbiguity {
			continue
		}
		if !m.IsCrosslinkRef() {
			defined[*label] = true
		}
	}
	for _, m := range mods {
		label := m.GetCrosslinkID()
		if label == nil || classifyLabel(*label) != labelKindAmbiguity {
			continue
		}
		if m.IsCrosslinkRef() && !defined[*label] {
			return fmt.Errorf("ambiguity group reference #%s has no matching definition", *label)
		}
	}
	return nil
}

// validateCrosslinkAndBranchLabels checks every crosslink/branch reference resolves to a
// definition. A definition with no reference is valid (a "dead end" crosslink).
func validateCrosslinkAndBranchLabels(mods []*Modification) error {
	definedCrosslinks := make(map[string]bool)
	hasBranchDefinition := false
	for _, m := range mods {
		if label := m.GetCrosslinkID(); label != nil && classifyLabel(*label) == labelKindCrosslink && !m.IsCrosslinkRef() {
			definedCrosslinks[*label] = true
		}
		if m.HasBranch() && !m.IsBranchRef() {
			hasBranchDefinition = true
		}
	}

	for _, m := range mods {
		if label := m.GetCrosslinkID(); label != nil && classifyLabel(*label) == labelKindCrosslink && m.IsCrosslinkRef() {
			if !definedCrosslinks[*label] {
				return fmt.Errorf("crosslink reference #%s has no matching definition", *label)
			}
		}
		if m.IsBranchRef() && !hasBranchDefinition {
			return fmt.Errorf("branch reference #BRANCH has no matching definition")
		}
	}
	return nil
}

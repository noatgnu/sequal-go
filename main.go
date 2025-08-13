package main

import (
	"fmt"
	"github.com/noatgnu/sequal-go/sequal"
)

func main() {
	fmt.Println("Sequal - Go ProForma Library")
	fmt.Println("============================")

	// Example usage of the library
	examples := []string{
		"PEPTIDE",
		"ELVIS[Phospho]K",
		"[Acetyl]-PEPTIDE-[Amidated]",
		"PEPTIDE/2",
		"PEPTIDE/2[+Na+]",
		"<[Carbamidomethyl]@C>PEPTCDE",
		"PEPTIDE/2+ANOTHER/3", // Chimeric
		"ELVIS[U:Phospho|+79.966331]K", // Complex modification
	}

	fmt.Println("\nParsing ProForma examples:")
	fmt.Println("--------------------------")

	for _, proforma := range examples {
		fmt.Printf("\nProForma: %s\n", proforma)
		
		seq, err := sequal.FromProforma(proforma)
		if err != nil {
			fmt.Printf("Error: %v\n", err)
			continue
		}

		fmt.Printf("Stripped sequence: %s\n", seq.ToStrippedString())
		fmt.Printf("Length: %d\n", seq.GetLength())
		
		if charge := seq.GetCharge(); charge != nil {
			fmt.Printf("Charge: %d\n", *charge)
		}
		
		if species := seq.GetIonicSpecies(); species != nil {
			fmt.Printf("Ionic species: %s\n", *species)
		}
		
		if seq.IsChimeric() {
			fmt.Printf("Chimeric peptidoforms: %d\n", len(seq.GetPeptidoforms()))
		}
		
		if seq.IsMultiChain() {
			fmt.Printf("Multi-chain sequences: %d\n", len(seq.GetChains()))
		}
		
		// Show global modifications
		globalMods := seq.GetGlobalMods()
		if len(globalMods) > 0 {
			fmt.Printf("Global modifications: %d\n", len(globalMods))
		}
		
		// Regenerate ProForma
		regenerated := seq.ToProforma()
		fmt.Printf("Regenerated: %s\n", regenerated)
	}

	fmt.Println("\n\nLibrary features demonstrated:")
	fmt.Println("- ProForma parsing and generation")
	fmt.Println("- Terminal modifications")
	fmt.Println("- Charge state handling")
	fmt.Println("- Global modifications")
	fmt.Println("- Chimeric sequences")
	fmt.Println("- Complex modification notation")
}
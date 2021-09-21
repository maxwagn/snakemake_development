from PyPDF2 import PdfFileReader, PdfFileWriter

CHROMOSOMES = [str(i) for i in range(1,25) if i not in [23] ] # +  ['MT'] # Mitogenome is nested within scaffold_245_arrow_ctg1


def merge_pdfs(paths, output):
    pdf_writer = PdfFileWriter()

    for path in paths:
        pdf_reader = PdfFileReader(path)
        for page in range(pdf_reader.getNumPages()):
            # Add each page to the writer object
            pdf_writer.addPage(pdf_reader.getPage(page))

    # Write out the merged PDF
    with open(output, 'wb') as out:
        pdf_writer.write(out)

if __name__ == '__main__':
    files = []
    for chr in CHROMOSOMES:
        files.append("202007Gouania.site_stats.chr{}.pdf".format(chr))
    merge_pdfs(files, output='Site_Stats_AllChrom.pdf')

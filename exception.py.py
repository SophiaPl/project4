from Bio.PDB import *
from Bio.SeqUtils.ProtParam import ProteinAnalysis, ProtParamData

try:
    def protein_seq_creator(pdb_name: str):
        try:
            pdbl = PDBList()
        except Exception:
            return 'Добавляйте библиотеку!'
        try:
            pdbl.retrieve_pdb_file(pdb_name)  # скачивание структуры белка в формате PDBx/mmCif
        except Exception:
            return 'Файл в другом формате!'

        try:
            path = str(pdb_name[1:3].lower() + '/' + pdb_name.lower() + '.cif')  # путь до скаченного файла
        except Exception:
            return 'Не правильно путь до файла!'
        parser = MMCIFParser()
        structure = parser.get_structure(pdb_name, path)  # чтение структуры из файла

        ppb = PPBuilder()  # чтение последовательности белка
        seq = ''
        for pp in ppb.build_peptides(structure):
            seq += pp.get_sequence()

        return seq
except Exception:
    print('Введен много аргументов!')
else:
    print('В первом функции не существует исключение!')

try:
    def surface_analysis(sequence: str):  # анализ структуры белка
        prot = ProteinAnalysis(str(sequence))
        kd = prot.protein_scale(window=7, param_dict=ProtParamData.kd)  # Kyte & Doolittle index of hydrophobicity
        em = prot.protein_scale(window=7, param_dict=ProtParamData.em)  # Surface accessibility
        seq_cut = sequence[3:len(sequence) - 3]
        # парметр window=7 задаёт длину окружения, в пределах которого считается гидрофобность/доступность
        # поэтому первые и последние три остатка не анализируются
        try:
            seq_first3 = sequence[:3].upper()
            seq_last3 = sequence[-3:].upper()
        except Exception:
            return 'Индекс не правильно!'

        residual_parameters = []
        try:
            for a, h, s in zip(seq_cut, kd, em):
                residual_parameters.append((a, round(h, 3), round(s, 3)))  # остаток + его kd и em
        except Exception:
            return 'Исклчения в for!'

        residuals_availability = ''
        for i in residual_parameters:
            if i[1] <= 0 and i[2] >= 1:
                residuals_availability += i[0].upper()  # большими буквами отмечены остатки на поверхности белка
            else:
                residuals_availability += i[0].lower()  # маленькими - остатки внутри белка

        residuals_availability = seq_first3 + residuals_availability + seq_last3
        # возвращение первых и последних трёх остатков

        return residuals_availability
except Exception:
    print('Введен много аргументов!')
else:
    print('В втором функции не существует исключение!')

try:
    def site_seq_creator(prot_site: str):
        try:
            cut_index = prot_site.index('/')
        except NameError:
            return 'имя не существует!'
        site_seq_left, site_seq_right = prot_site[:cut_index], prot_site[cut_index + 1:]
        site_seq = (site_seq_left + site_seq_right).upper()
        return site_seq, site_seq_left, site_seq_right
        # последовательность сайта без /, часть последовательности до /, часть после /
except Exception:
    print('Введен много аргументов!')
else:
    print('В третьем функции не существует исключение!')



def proteolysis_sites(seq: str, site_seq: str, site_seq_left: str, site_seq_right: str):
    step = 0
    step_and_index = {}
    proteins = [seq]

    def proteolysis(proteins: list):
        for i in range(len(proteins)):
            proteins[i] = surface_analysis(proteins[i])
            # анализ всего белка и впоследствии его частей функцией surface_analysis()
        indexes = []  # координаты сайтов протеолиза
        new_proteins = []  # новые части белка после протеолиза
        for prot in proteins:
            if site_seq in prot:  # если сайт протеолиза есть в последовательности белка (большими буквами)
                id = 0
                for i in range(prot.count(site_seq)):
                    id += prot.index(site_seq) + len(site_seq_left)  # координата текущего сайта протеолиза
                    indexes.append(str(id) + '/' + str(id + 1))
                    new_part = (prot[:prot.index(site_seq)]).upper() + site_seq_left
                    # новая часть белка после протеолиза (всё, что левее /)
                    new_proteins.append(new_part)
                    prot = prot[prot.index(site_seq):]
                    prot = prot[prot.index(site_seq_right):]
                    # новая часть белка после протеолиза (всё, что правее /)
            else:
                continue
        return new_proteins, indexes

    while any([site_seq in prot for prot in proteins]):
        # пока есть хоть одна часть белка с последовательностью сайта протеолиза, работает функция proteolysis()
        proteins, indexes = proteolysis(proteins)
        step += 1
        step_and_index[step] = indexes  # шаг и найдеденные сайты протеолиза
    else:
        if step_and_index:
            return step_and_index
        else:
            return 'No sites of proteolysis were found'
            # если в белке в принципе нет последовательности сайта протеолиза


# временный блок, пока не пропишем cmd_args
try:
    if __name__ == '__main__':
        pdb = input('Input protein structure name on RCSB PDB: ')
        protease_site = input('Input protease cleavage site in AA/AA format: ')

        pseq = protein_seq_creator(pdb)
        sseq, sseql, sseqr = site_seq_creator(protease_site)

        print(proteolysis_sites(pseq, sseq, sseql, sseqr))
except FileNotFoundError:
    print('Не существует файла в базза данных!')
finally:
    print("End")
% Прочитать файл
function input_bits = file_reader(File, Size_Buffer)
    % Прочитать черно-белую картинку
    grayImage=imread(File);
    
    % Перевести картинку в битовый формат
    binaryImage = imbinarize(grayImage);
    input_bits = double(binaryImage); % из логического массива в массив 0 и 1

    % Взять нужное количество бит из всей картинки
    input_bits = input_bits(1:Size_Buffer);
    
end
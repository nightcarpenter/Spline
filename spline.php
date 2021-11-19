<?php


$width = 850;
$height = 850;

// Можно задать толщину линии и точек
$thickness = 1;
$radius = 6;

// Можно задать произвольное количество точек
$x = array(50, 175, 220, 329, 450, 531, 625, 710, 799);
$y = array(200, 450, 120, 587, 650, 415, 567,79, 699);
$n = count($x);

$img = imageCreateTrueColor($width, $height);
$background_color = imageColorAllocate($img, 235, 235, 235);
$black = imageColorAllocate($img, 0, 0, 0);
$grey = imageColorAllocate($img, 0, 0, 255);
imageFilledRectangle($img, 0, 0, $width, $height, $background_color);




// Отметим точки на чертеже
for($i=0; $i<count($x); $i++){
	imageFilledEllipse($img, X($x[$i]), Y($y[$i]), $radius, $radius, $black);
}

// Составим систему уравнений для коэффициентов c
$SLAE = CreateMatrixC($x, $y, $n);

// Решаем систему уравнений
$solution = SolvingSLAE($SLAE);

// Добавляем первый нулевой элемент (c0 = 0)
$c = array();
$c[0] = 0;
for($i = 0; $i < count($solution); $i++){
	$c[$i+1] = $solution[$i];
}
$c[count($solution)+1] = 0;

// Создадим массив из коэффициентов для кубических многочленов
$ABCD = CreateMatrixABCD($x, $y, $n, $c);

// Рисуем сплайн
for ($i = 0; $i < count($ABCD); $i++){
	$f = $ABCD[$i];
	DrawSpline($f, $x[$i], $x[$i+1]);
}




imageJpeg($img, "spline.jpg");
imagedestroy($img);
echo '<img src="spline.jpg">';




// Составим систему уравнений для коэффициентов c
function CreateMatrixC($x, $y, $n){
	// Создаём матрицу из нулей нужного размера
	$SLAE = array();
	for($j = 0; $j < ($n - 2); $j++){
		$SLAE[$j] = array();
		for($i = 0; $i <= $n; $i++){
			$SLAE[$j][$i] = 0;
		}
	}
	
	// Создаём массив h
	$h = CreateH($x, $n);
	
	//Заполняем матрицу коэффициентами
	for($j = 0; $j < ($n - 2); $j++){
		if($j == 0){
			$i = 1;
			$SLAE[$j][$i] = 2/3*($h[$i-1] + $h[$i]);
			$SLAE[$j][$i+1] = $h[$i]/3;
			$SLAE[$j][$n] = ($y[$i+1] - $y[$i])/$h[$i] - ($y[$i] - $y[$i-1])/$h[$i-1];
		}
		else if ($j == ($n - 3)){
			$i = $n - 2;
			$SLAE[$j][$i-1] = $h[$i-1]/3;
			$SLAE[$j][$i] = 2/3*($h[$i-1] + $h[$i]);
			$SLAE[$j][$n] = ($y[$i+1] - $y[$i])/$h[$i] - ($y[$i] - $y[$i-1])/$h[$i-1];
		}
		else{
			$i = $j + 1;
			$SLAE[$j][$i-1] = $h[$i-1]/3;
			$SLAE[$j][$i] = 2/3*($h[$i-1] + $h[$i]);
			$SLAE[$j][$i+1] = $h[$i]/3;
			$SLAE[$j][$n] = ($y[$i+1] - $y[$i])/$h[$i] - ($y[$i] - $y[$i-1])/$h[$i-1];
		}
	}
	
	//Убираем первую и последнюю нулевые строки
	$SLAE_MOD = array();
	for($j = 0; $j < count($SLAE); $j++){
		$SLAE_MOD[$j] = array();
		$k = 0;
		for($i = 0; $i < count($SLAE[$j]); $i++){
			if($i != 0 and $i != count($SLAE[$j]) - 2){
				$SLAE_MOD[$j][$k] = $SLAE[$j][$i];
				$k += 1;
			}
		}
	}
	
	return $SLAE_MOD;
}




// Создаём массив h
function CreateH($x, $n){
	$h = array();
	for($i = 0; $i < $n - 1; $i++){
		$h[$i] = $x[$i+1] - $x[$i];
	}
	return $h;
}




// Для решения СЛАУ
function SolvingSLAE($M){
	$nj = count($M);
	$ni = count($M[0]);
	
	$y = array();
	$x = array();
	
	$a = array();
	$b = array();
	
	// Прямой ход
	$y[0] = $M[0][0];
	$a[0] = -$M[0][1]/$y[0];
	$b[0] = $M[0][$ni-1]/$y[0];
	
	for($i = 1; $i < $nj; $i++){
		$y[$i] = $M[$i][$i] + ($M[$i][$i-1] * $a[$i-1]);
		$a[$i] = -$M[$i][$i+1]/$y[$i];
		$b[$i] = ($M[$i][$ni-1] - $M[$i][$i-1]*$b[$i-1])/$y[$i];
	}
	
	// Обратный ход
	$x[$nj-1] = $b[$nj-1];
	
	for($i = ($nj - 2); $i >= 0; $i--){
		$x[$i] = $a[$i]*$x[$i+1] + $b[$i];
	}
	
	return $x;
}




// Для создания матрицы всех коэффициентов кубических многочленов
function CreateMatrixABCD($x, $y, $n, $c){
	$h = CreateH($x, $n);
	
	// Массив коэффициентов a
	$a = array();
	for($i = 0; $i < $n; $i++){
		$a[$i] = $y[$i];
	}
	
	// Массив коэффициентов b и d
	$d = array();
	$b = array();
	for($i = 0; $i < $n - 1; $i++){
		$d[$i] = ($c[$i+1] - $c[$i])/(3*$h[$i]);
		
		$b[$i] = ($y[$i+1] - $y[$i])/$h[$i] - ($h[$i]/3)*($c[$i+1] + 2*$c[$i]);
	}
	
	// Создаём матрицу со всеми коэффициентами
	$ABCD = array();
	
	for($i = 0; $i < count($c)-1; $i++){
		$ABCD[$i][0] = $a[$i];
		$ABCD[$i][1] = $b[$i];
		$ABCD[$i][2] = $c[$i];
		$ABCD[$i][3] = $d[$i];
	}
	
	return $ABCD;
}




// Построения сплайна между 3-мя точками
function DrawSpline($f, $x_min, $x_max){
	global $img, $black;
	
	$s = ($x_max - $x_min) / 200;
	$x1 = $x_min;
	
	do{
		$x2 = $x1 + $s;
		
		$y1 = F($f, $x1, $x_min);
		$y2 = F($f, $x2, $x_min);
		
		imageLine($img, X($x1), Y($y1), X($x2), Y($y2), $black); 
		
		$x1 = $x2;
	}while($x1 < $x_max);
}


// Функция кубического четырёхчлена для построения сплайна
function F($f, $x, $x0){
	$y = $f[0] + $f[1]*($x-$x0) + $f[2]*pow(($x-$x0), 2) + $f[3]*pow(($x-$x0), 3);
	return $y;
}




// Преобразование координат к пиксельным
function X($x){
	return $x;
}

function Y($y){
	global $height;
	return $height - $y;
}


// Вывод матрицы на экран (для проверки)
function ShowMatrix($M){
	echo '</br>';
	for($j = 0; $j < count($M); $j++){
		for($i = 0; $i < count($M[$j]); $i++){
			echo round($M[$j][$i], 100), ' / ';
		}
		echo '</br>';
	}
	echo '</br>';
}


// Вывод массива на экран (для проверки)
function ShowArray($arr){
	echo '</br>';
	for($i = 0; $i < count($arr); $i++){
		echo $arr[$i];
		if($i != count($arr) - 1){
			echo ' / ';
		}
	}
	echo '</br>';
}


?>
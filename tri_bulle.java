public class InsertionSort {

    public static void insertionSort(int[] array) {
        for (int i = 1; i < array.length; i++) {
            int key = array[i];      // Élément à insérer
            int j = i - 1;

            // Décale les éléments plus grands que key
            while (j >= 0 && array[j] > key) {
                array[j + 1] = array[j];
                j--;
            }

            // Insère key à la bonne position
            array[j + 1] = key;
        }
    }

    public static void main(String[] args) {
        int[] numbers = {5, 2, 9, 1, 5, 6};

        insertionSort(numbers);

        System.out.println("Tableau trié :");
        for (int num : numbers) {
            System.out.print(num + " ");
        }
    }
}
